//
//  Advection.ispc
//  CubismUP_3D
//
//  Created by Christian Conti on 2/4/16.
//  Copyright © 2016 ETHZ. All rights reserved.
//

#ifndef _SP_COMP_
typedef double Real;
#else // _SP_COMP_
typedef float Real;
#endif // _SP_COMP_

// Euler only
export void advect_ispc(uniform Real uOut[],
						uniform Real vOut[],
						uniform Real wOut[],
						uniform Real uOutRHS[], //tmp
						uniform Real vOutRHS[],
						uniform Real wOutRHS[],
						const uniform Real uIn[], //lab
						const uniform Real vIn[],
						const uniform Real wIn[],
#ifdef _MOVING_FRAME_
						uniform Real uBody,
						uniform Real vBody,
						uniform Real wBody,
#endif
						const uniform Real factor)
{
	enum {
		_BSX2_ = _BSX_+4,
		_BSY2_ = _BSY_+4,
		_BSZ2_ = _BSZ_+4,
		SLICESIZE = _BSX_*_BSY_,
		SLICESIZE2 = (_BSX_+4)*(_BSY_+4)
	};
	
	cfor (uniform int k=0; k<_BSZ_; k++)
	{
		uniform int baseZ = k*SLICESIZE2;
		uniform int baseZout = k*SLICESIZE;
		
		cfor (uniform int j=0; j<_BSY_; j++)
		{
			uniform int base = baseZ + j*_BSX2_;
			uniform int baseOut = baseZout + j*_BSX_;
			
			// 150 FLOP
			foreach (i=0 ... _BSX_)
			{
				// 3 Reads
#ifndef _MOVING_FRAME_
				const Real u = uOut[baseOut+i];
				const Real v = vOut[baseOut+i];
				const Real w = wOut[baseOut+i];
#else
				// 3 FLOP
				const Real u = uOut[baseOut+i] - uBody;
				const Real v = vOut[baseOut+i] - vBody;
				const Real w = wOut[baseOut+i] - wBody;
#endif
				
				// 6 FLOP
				const Real minu = min(u,(Real)0.);
				const Real maxu = max(u,(Real)0.);
				const Real minv = min(v,(Real)0.);
				const Real maxv = max(v,(Real)0.);
				const Real minw = min(w,(Real)0.);
				const Real maxw = max(w,(Real)0.);
				
				// u - 47 FLOP, 15 R/W
				// 1 FLOP, 1 Read
				const Real u3 = 3 * uIn[base + i];
				
				// 11 FLOP, 4 Read
				const Real uW  = uIn[base + i-1];
				const Real uE  = uIn[base + i+1];
				const Real uW2 = uIn[base + i-2];
				const Real uE2 = uIn[base + i+2];
				const Real dudxW =  2*uE  + u3 - 6*uW +   uW2;
				const Real dudxE = -  uE2 + 6*uE - u3 - 2*uW;
				
				// 11 FLOP, 4 Read
				const Real uS  = uIn[base + i-_BSX2_];
				const Real uN  = uIn[base + i+_BSX2_];
				const Real uS2 = uIn[base + i-2*_BSX2_];
				const Real uN2 = uIn[base + i+2*_BSX2_];
				const Real dudyS = 2*uN  + u3 - 6*uS +   uS2;
				const Real dudyN = - uN2 + 6*uN - u3 - 2*uS;
				
				// 11 FLOP, 4 Read
				const Real uF  = uIn[base + i-SLICESIZE2];
				const Real uB  = uIn[base + i+SLICESIZE2];
				const Real uF2 = uIn[base + i-2*SLICESIZE2];
				const Real uB2 = uIn[base + i+2*SLICESIZE2];
				const Real dudzF = 2*uB  + u3 - 6*uF +   uF2;
				const Real dudzB = - uB2 + 6*uB - u3 - 2*uF;
				
				// 13 FLOP, 1 Read, 1 Write
				uOutRHS[baseOut+i] = uOut[baseOut+i] + factor*(maxu * dudxW + minu * dudxE +
															   maxv * dudyS + minv * dudyN +
															   maxw * dudzF + minw * dudzB);
				
				// v
				const Real v3 = 3 * vIn[base + i];
				
				const Real vW  = vIn[base + i-1];
				const Real vE  = vIn[base + i+1];
				const Real vW2 = vIn[base + i-2];
				const Real vE2 = vIn[base + i+2];
				const Real dvdxW =  2*vE  + v3 - 6*vW +   vW2;
				const Real dvdxE = -  vE2 + 6*vE - v3 - 2*vW;
				
				const Real vS  = vIn[base + i-_BSX2_];
				const Real vN  = vIn[base + i+_BSX2_];
				const Real vS2 = vIn[base + i-2*_BSX2_];
				const Real vN2 = vIn[base + i+2*_BSX2_];
				const Real dvdyS = 2*vN  + v3 - 6*vS +   vS2;
				const Real dvdyN = - vN2 + 6*vN - v3 - 2*vS;
				
				const Real vF  = vIn[base + i-SLICESIZE2];
				const Real vB  = vIn[base + i+SLICESIZE2];
				const Real vF2 = vIn[base + i-2*SLICESIZE2];
				const Real vB2 = vIn[base + i+2*SLICESIZE2];
				const Real dvdzF = 2*vB  + v3 - 6*vF +   vF2;
				const Real dvdzB = - vB2 + 6*vB - v3 - 2*vF;
				
				vOutRHS[baseOut+i] = vOut[baseOut+i] + factor*(maxu * dvdxW + minu * dvdxE +
															   maxv * dvdyS + minv * dvdyN +
															   maxw * dvdzF + minw * dvdzB);
				
				// w
				const Real w3 = 3 * wIn[base + i];
				
				const Real wW  = wIn[base + i-1];
				const Real wE  = wIn[base + i+1];
				const Real wW2 = wIn[base + i-2];
				const Real wE2 = wIn[base + i+2];
				const Real dwdxW =  2*wE  + w3 - 6*wW +   wW2;
				const Real dwdxE = -  wE2 + 6*wE - w3 - 2*wW;
				
				const Real wS  = wIn[base + i-_BSX2_];
				const Real wN  = wIn[base + i+_BSX2_];
				const Real wS2 = wIn[base + i-2*_BSX2_];
				const Real wN2 = wIn[base + i+2*_BSX2_];
				const Real dwdyS = 2*wN  + w3 - 6*wS +   wS2;
				const Real dwdyN = - wN2 + 6*wN - w3 - 2*wS;
				
				const Real wF  = wIn[base + i-SLICESIZE2];
				const Real wB  = wIn[base + i+SLICESIZE2];
				const Real wF2 = wIn[base + i-2*SLICESIZE2];
				const Real wB2 = wIn[base + i+2*SLICESIZE2];
				const Real dwdzF = 2*wB  + w3 - 6*wF +   wF2;
				const Real dwdzB = - wB2 + 6*wB - w3 - 2*wF;
				
				wOutRHS[baseOut+i] = wOut[baseOut+i] + factor*(maxu * dwdxW + minu * dwdxE +
															   maxv * dwdyS + minv * dwdyN +
															   maxw * dwdzF + minw * dwdzB);
			}
		}
	}
}