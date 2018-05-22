#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;

namespace FakeShader {
	
	float saturate(float x){ return std::max(0.f, std::min(1.f, x)); }
	
	float smoothstep(float edge0, float edge1, float x)
	{
		// Scale, bias and saturate x to 0..1 range
		x = saturate((x - edge0)/(edge1 - edge0)); 
		// Evaluate polynomial
		return x*x*(3 - 2*x);
	}
	
	// red
	const float cr = .45;
	const float wr = .075;
	
	// green
	const float cg = .6;
	const float wg = .05;
	
	// blue
	const float cb = .2;//.15;
	const float wb = .1;//.075;

	// affects alpha
	float modulate_material(float rho)
	{	
		float boxr = (1.-smoothstep(cr, cr+wr, rho))*smoothstep(cr-wr, cr, rho);
		float boxg = (1.-smoothstep(cg, cg+wg, rho))*smoothstep(cg-wg, cg, rho);
		float boxb = (1.-smoothstep(cb, cb+wb, rho))*smoothstep(cb-wb, cb, rho);
		return 1*boxb + 2*boxr + 8*boxg;
		//5.*boxr + 2.5*boxg + 1.2*boxb;
		//return 1.*boxr + 1.*boxg + 1.*boxb;
	}
	
	// affects gaussian-box (small (<1): gaussian, large: box)
	float modulate_lighting_effect_red(float rho)
	{
		float boxr = (1.-smoothstep(cr, cr+wr, rho))*smoothstep(cr-wr, cr, rho);
		return .5*boxr;
	}
	
	float modulate_lighting_effect_green(float rho)
	{
		float boxg = (1.-smoothstep(cg, cg+wg, rho))*smoothstep(cg-wg, cg, rho);
		return 1.*boxg;
	}
	
	float modulate_lighting_effect_blue(float rho)
	{
		float boxb = (1.-smoothstep(cb, cb+wb, rho))*smoothstep(cb-wb, cb, rho);
		return 1*boxb;//.3*boxb;
	}
	
	typedef float RGBA[4];
	
	struct vec3
	{
		float data[3];
		
		vec3(const float x)
		{
			for(int i = 0; i< 3; ++i)
				data[i] = x;
		}
		
		vec3(const float x, const float y, const float z)
		{
			data[0] = x;
			data[1] = y;
			data[2] = z;
		}
	};
	
	vec3 normalize(vec3 x)
	{
		const double IxI2 = x.data[0]*x.data[0] + x.data[1]*x.data[1] + x.data[2]*x.data[2];
		for(int i=0; i<3; ++i)
			x.data[i] /= IxI2 + 1e-6;
		
		return x;
	}
	
	vec3 min(vec3 x, vec3 y)
	{
		for(int i = 0; i < 3; ++i)
			x.data[i] = std::min(x.data[i], y.data[i]);
		
		return x;
	}
	
	vec3 operator *(vec3 x, vec3 y)
	{
		for(int i = 0; i < 3; ++i)
			x.data[i] *= y.data[i];
		
		return x;
	}

	vec3 operator +(vec3 x, vec3 y)
	{
		for(int i = 0; i < 3; ++i)
			x.data[i] += y.data[i];
		
		return x;
	}
	
	vec3 operator *(vec3 x, float a)
	{
		for(int i = 0; i < 3; ++i)
			x.data[i] *= a;
		
		return x;
	}
	
	vec3 operator *(float a, vec3 x)
	{
		for(int i = 0; i < 3; ++i)
			x.data[i] *= a;
		
		return x;
	}
	
	void shader_main(float rho, RGBA& color)
	{
		const float factor = 5;//15.0;
		color[0] = 0;
		color[1] = 0;
		color[2] = 0;
		color[3] = std::min(1.f, modulate_material(rho))*factor;
		
		// min(vec3(maxvalue), normalize(color)*modulate) + min(vec3(d * modulate), vec3(e + f * pow(max(0.f,e),g)));
		// maxvalue: maximum value that the color can assume, determines maximum intensity
		// color: base color assigned to the given gaussian
		// modulate: modulation of the color - gaussian/box
		// d: (specular)
		// cos_value: set to 0, manages specular effects (specular)
		// f: (specular)
		// g: (specular)
		
		// red
		const float cos_value = 0;
		float redmodulate = modulate_lighting_effect_red(rho);		
		vec3 cazzo = min(vec3(.4), normalize(vec3(.1,0,.0))*redmodulate) + min(vec3(0.2*redmodulate), vec3(cos_value + 0.5*pow(max(0.f,cos_value),10.f)));
		for(int i=0; i<3; ++i) color[i] += cazzo.data[i]*factor;
		
		// green
		float greenmodulate = modulate_lighting_effect_green(rho);
		cazzo = min(vec3(.3), normalize(vec3(1.0,1.0,0.0))*greenmodulate) + min(vec3(.01*greenmodulate), vec3(cos_value + 0.5*pow(max(0.f,cos_value),10.f)));
		for(int i=0; i<3; ++i) color[i] += 10*cazzo.data[i]*factor;
		
		// blue - purple
		float bluemodulate = modulate_lighting_effect_blue(rho);
		cazzo = min(vec3(.6), normalize(vec3(.0,.15,1.))*bluemodulate) + min(vec3(1.*bluemodulate), vec3(cos_value + 0.5*pow(max(0.f,cos_value),10.f)));
		for(int i=0; i<3; ++i) color[i] += cazzo.data[i]*factor;
		
		// saturate/desaturate
		for (int i = 0; i < 3; ++i)
			color[i] *= 2.5;
		
	}

	vector<float> generate_lut(const int N)
	{
		vector<float> lut(4*N);
		
		const double h = 1./N;
		
		for(int i = 0; i < N ; ++i)
		{
			const float rho = (0.5 + i) * h;
			
			RGBA color;
			shader_main(rho, color);
			
			lut[4*i + 0] = color[0];
			lut[4*i + 1] = color[1];
			lut[4*i + 2] = color[2];
			lut[4*i + 3] = color[3];
		}
		
		return lut;
	}
}

int main (int argc, char * const argv[]) {
	
	const int N =  4096;
	vector<float> mylut = FakeShader::generate_lut(N);
	
	FILE * f = fopen("mylut.txt", "w");
	for(int i = 0 ; i < N; ++i)
		fprintf(f, "%i %e %e %e %e\n", i, mylut[4*i], mylut[4*i + 1], mylut[4*i + 2], mylut[4*i + 3]);
	
	fclose(f);
	
    // insert code here...
    std::cout << "Hello, World!\n";
    return 0;
}

