//
//  IF3D_CarlingFishOperator.cpp
//  IncompressibleFluids3D
//
//  Created by Wim van Rees on 4/15/13.
//
//

#include "IF3D_CarlingFishOperator.h"
#include "GenericOperator.h"

struct VolumeSegment_OBB
{
    std::pair<int, int> s_range;
    Real normalI[3]; // should be normalized and >=0
    Real normalJ[3];
    Real normalK[3];
    Real w[3]; // halfwidth
    Real c[3]; // center

#if 1==0
    VolumeSegment_OBB(std::pair<int, int> s_range, const Real bbox[3][2])
    : s_range(s_range)
    {
        normalI[1]=normalI[2]=normalJ[0]=normalJ[2]=normalK[0]=normalK[1]=0.0;
        normalI[0]=normalJ[1]=normalK[2]=1.0;
        for(int i=0;i<3;++i) {
        	w[i] = 0.5*(bbox[i][1]-bbox[i][0]);
        	c[i] = bbox[i][0] + w[i];
        	assert(w[i]>0);
        }

    }
#endif
    
    VolumeSegment_OBB() { }

    void prepare(std::pair<int, int> _s_range, const Real bbox[3][2])
	{
    	normalI[1]=normalI[2]=normalJ[0]=normalJ[2]=normalK[0]=normalK[1]=0.0;
		normalI[0]=normalJ[1]=normalK[2]=1.0;
    	s_range.first = _s_range.first;
    	s_range.second = _s_range.second;
		for(int i=0;i<3;++i) {
			w[i] = 0.5*(bbox[i][1]-bbox[i][0]);
			c[i] = bbox[i][0] + w[i];
			assert(w[i]>0);
		}
	}

    void normalizeNormals()
    {
        const Real magI = std::sqrt(normalI[0]*normalI[0]+normalI[1]*normalI[1]+normalI[2]*normalI[2]);
        const Real magJ = std::sqrt(normalJ[0]*normalJ[0]+normalJ[1]*normalJ[1]+normalJ[2]*normalJ[2]);
        const Real magK = std::sqrt(normalK[0]*normalK[0]+normalK[1]*normalK[1]+normalK[2]*normalK[2]);
        assert(magI > std::numeric_limits<Real>::epsilon());
        assert(magJ > std::numeric_limits<Real>::epsilon());
        assert(magK > std::numeric_limits<Real>::epsilon());
        const Real invMagI = 1.0/magI;
        const Real invMagJ = 1.0/magJ;
        const Real invMagK = 1.0/magK;
        
        for(int i=0;i<3;++i) {
	    // also take absolute value since thats what we need when doing intersection checks later
            normalI[i]=std::abs(normalI[i])*invMagI;
            normalJ[i]=std::abs(normalJ[i])*invMagJ;
            normalK[i]=std::abs(normalK[i])*invMagK;
        }
    }

    void changeToComputationalFrame(const Real position[3], const Real quaternion[4])
    {
        // we are in CoM frame and change to comp frame --> first rotate around CoM (which is at (0,0) in CoM frame), then update center
        const Real w = quaternion[0];
        const Real x = quaternion[1];
        const Real y = quaternion[2];
        const Real z = quaternion[3];
        const Real Rmatrix3D[3][3] = {
            {1.-2*(y*y+z*z),  2*(x*y-z*w),    2*(x*z+y*w)},
            {2*(x*y+z*w),    1.-2*(x*x+z*z),  2*(y*z-x*w)},
            {2*(x*z-y*w),    2*(y*z+x*w),    1.-2*(x*x+y*y)}
        };
        const Real p[3] = {c[0],c[1],c[2]};
        const Real nx[3] = {normalI[0],normalI[1],normalI[2]};
        const Real ny[3] = {normalJ[0],normalJ[1],normalJ[2]};
        const Real nz[3] = {normalK[0],normalK[1],normalK[2]};
        for(int i=0;i<3;++i) {
            c[i] = Rmatrix3D[i][0]*p[0] + Rmatrix3D[i][1]*p[1] + Rmatrix3D[i][2]*p[2];
            normalI[i] = Rmatrix3D[i][0]*nx[0] + Rmatrix3D[i][1]*nx[1] + Rmatrix3D[i][2]*nx[2];
            normalJ[i] = Rmatrix3D[i][0]*ny[0] + Rmatrix3D[i][1]*ny[1] + Rmatrix3D[i][2]*ny[2];
            normalK[i] = Rmatrix3D[i][0]*nz[0] + Rmatrix3D[i][1]*nz[1] + Rmatrix3D[i][2]*nz[2];
        }
        c[0] +=position[0];
        c[1] +=position[1];
        c[2] +=position[2];
        normalizeNormals();
    }
    
    bool isIntersectingWithAABB(const Real start[3],const Real end[3], const Real safe_distance = 0.0) const
    {
    	//start and end are two diagonally opposed corners of grid block
        const Real AABB_w[3] = { //half block width + safe distance
            0.5*(end[0] - start[0]) + 2.0*safe_distance,
            0.5*(end[1] - start[1]) + 2.0*safe_distance,
            0.5*(end[2] - start[2]) + 2.0*safe_distance
        };
        const Real AABB_c[3] = { //block center
            start[0] + AABB_w[0] - safe_distance,
            start[1] + AABB_w[1] - safe_distance,
            start[2] + AABB_w[2] - safe_distance
        };
        assert(AABB_w[0]>0);
        assert(AABB_w[1]>0);
        assert(AABB_w[2]>0);
        bool intersects = true;
        Real r;
        {
            r = w[0]*normalI[0] + w[1]*normalJ[0] + w[2]*normalK[0];
            intersects &= ((c[0]-r <= AABB_c[0] + AABB_w[0]) && (c[0]+r >= AABB_c[0] - AABB_w[0]));
            
            r = w[0]*normalI[1] + w[1]*normalJ[1] + w[2]*normalK[1];
            intersects &= ((c[1]-r <= AABB_c[1] + AABB_w[1]) && (c[1]+r >= AABB_c[1] - AABB_w[1]));
            
            r = w[0]*normalI[2] + w[1]*normalJ[2] + w[2]*normalK[2];
            intersects &= ((c[2]-r <= AABB_c[2] + AABB_w[2]) && (c[2]+r >= AABB_c[2] - AABB_w[2]));
        }
        {
            //r = AABB_w[0]*normalI[0] + AABB_w[1]*normalI[1] + AABB_w[2]*normalI[2];
            r = AABB_w[0]*normalI[0] + AABB_w[1]*normalJ[0] + AABB_w[2]*normalK[0];
            intersects &= ((AABB_c[0]-r <= c[0] + w[0]) && (AABB_c[0]+r >= c[0] - w[0]));
            
            //r = AABB_w[0]*normalJ[0] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalJ[2];
            r = AABB_w[0]*normalI[1] + AABB_w[1]*normalJ[1] + AABB_w[2]*normalK[1];
            intersects &= ((AABB_c[1]-r <= c[1] + w[1]) && (AABB_c[1]+r >= c[1] - w[1]));
            
            //r = AABB_w[0]*normalK[0] + AABB_w[1]*normalK[1] + AABB_w[2]*normalK[2];
            r = AABB_w[0]*normalI[2] + AABB_w[1]*normalJ[2] + AABB_w[2]*normalK[2];
            intersects &= ((AABB_c[2]-r <= c[2] + w[2]) && (AABB_c[2]+r >= c[2] - w[2]));
        }
        return intersects;
    }
};

struct PutFishOnBlocks
{    
    const Fish::FishMidlineData * cfish;
    Real position[3];
    Real quaternion[4];
    Real Rmatrix3D[3][3];

    PutFishOnBlocks(const Fish::FishMidlineData* const cfish, const Real pos[3], const Real quat[4]):
	cfish(cfish)
	{
        position[0]=pos[0];
        position[1]=pos[1];
        position[2]=pos[2];
        quaternion[0]=quat[0];
        quaternion[1]=quat[1];
        quaternion[2]=quat[2];
        quaternion[3]=quat[3];
        computeRotationMatrix();
	}
    
    void computeRotationMatrix()
    {
        const Real w = quaternion[0];
        const Real x = quaternion[1];
        const Real y = quaternion[2];
        const Real z = quaternion[3];
        const Real R[3][3] = {
            {1-2*(y*y+z*z),  2*(x*y-z*w),    2*(x*z+y*w)},
            {2*(x*y+z*w),    1-2*(x*x+z*z),  2*(y*z-x*w)},
            {2*(x*z-y*w),    2*(y*z+x*w),    1-2*(x*x+y*y)}
        };
        memcpy(Rmatrix3D, R, sizeof(R));
    }
    
    inline int find_closest_dist(const int s, const int dir, const Real x[3], Real & oldDistSq) const
    {        
        if((s+dir)<cfish->iFishStart or (s+dir)>cfish->iFishEnd)
            return s;
        
        const Real newDistSq = (x[0]-cfish->rX[s+dir])*(x[0]-cfish->rX[s+dir])
        		               + (x[1]-cfish->rY[s+dir])*(x[1]-cfish->rY[s+dir])
							   + (x[2])*(x[2]);
        
        if(oldDistSq<=newDistSq) {
            return s;
        } else {
            oldDistSq = newDistSq;
            return s+dir;
        }
    }
    
    void changeVelocityToComputationalFrame(Real x[3]) const
    {
        const Real p[3] = {x[0],x[1],x[2]};
        
        // rotate (around CoM)
        x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[0][1]*p[1] + Rmatrix3D[0][2]*p[2];
        x[1]=Rmatrix3D[1][0]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[1][2]*p[2];
        x[2]=Rmatrix3D[2][0]*p[0] + Rmatrix3D[2][1]*p[1] + Rmatrix3D[2][2]*p[2];
    }
    
    template<typename T>
    void changeToComputationalFrame(T x[3]) const
    {
        const T p[3] = {x[0],x[1],x[2]};
        
        // rotate (around CoM)
        x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[0][1]*p[1] + Rmatrix3D[0][2]*p[2];
        x[1]=Rmatrix3D[1][0]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[1][2]*p[2];
        x[2]=Rmatrix3D[2][0]*p[0] + Rmatrix3D[2][1]*p[1] + Rmatrix3D[2][2]*p[2];
        
        // translate
        x[0]+=position[0];
        x[1]+=position[1];
        x[2]+=position[2];
    }
    
    template<typename T>
    void changeFromComputationalFrame(T x[3]) const
    {
        // translate back to CoM
        const T p[3] = {
            x[0]-position[0],
            x[1]-position[1],
            x[2]-position[2]
        };
        
        // rotate back around CoM
        x[0]=Rmatrix3D[0][0]*p[0] + Rmatrix3D[1][0]*p[1] + Rmatrix3D[2][0]*p[2];
        x[1]=Rmatrix3D[0][1]*p[0] + Rmatrix3D[1][1]*p[1] + Rmatrix3D[2][1]*p[2];
        x[2]=Rmatrix3D[0][2]*p[0] + Rmatrix3D[1][2]*p[1] + Rmatrix3D[2][2]*p[2];
    }

    Real getSmallerDistToMidline(const int start_s, const Real x[3], int & final_s) const
    {
        Real relX[3] = {x[0],x[1],x[2]};
        changeFromComputationalFrame(relX);
        
        const Real curDistSq =  std::pow(relX[0]-cfish->rX[start_s],2)
        					    + std::pow(relX[1]-cfish->rY[start_s],2)
        						+ std::pow(relX[2],2);
        
        Real distSq;

        distSq = curDistSq; // check right
        const int sRight = find_closest_dist(start_s, +1, relX, distSq);

        distSq = curDistSq; // check left
        const int sLeft = find_closest_dist(start_s, -1, relX, distSq);

        if(sRight==start_s and sLeft==start_s) {
            final_s = start_s;
            return distSq;
        }
        
        assert(sRight==start_s or sLeft==start_s);
        
        int curr_s = start_s;
        int new_s = sRight == start_s ? sLeft : sRight;
        const int dir = new_s-curr_s;
        while(curr_s not_eq new_s) {
            curr_s = new_s;
            new_s = find_closest_dist(curr_s,dir,relX, distSq);
        }

        final_s = new_s;
        return distSq;
    }

	//inline void operator()(const BlockInfo& info, FluidBlock3D& b) const
	void operator()(const BlockInfo& info, FluidBlock& b, ObstacleBlock* const defblock, const std::vector<VolumeSegment_OBB>& vSegments) const
	{
        Real org[3];
        info.pos(org, 0, 0, 0);        
        const Real invh = 1.0/info.h_gridpoint;
        
        {
            Real * const chi = &(defblock->chi[0][0][0]);
            Real * const udef = &(defblock->udef[0][0][0][0]);
            
            static const int n = FluidBlock::sizeZ*FluidBlock::sizeY*FluidBlock::sizeX;
            for(int i=0; i<n; i++) {
                chi[i]=-1;
                udef[3*i+0]=0;
                udef[3*i+1]=0;
                udef[3*i+2]=0;
            }
        }
        
        // construct the shape (P2M with min(distance) as kernel) onto defblocks
        for(int i=0;i<vSegments.size();++i) {
        	//iterate over segments contained in the vSegm intersecting this block:
        	const int firstSegm = std::max(vSegments[i].s_range.first,cfish->iFishStart);
        	const int lastSegm = std::min(vSegments[i].s_range.second, cfish->iFishEnd);
        	for(int ss=firstSegm; ss<=lastSegm; ++ss) {
        		assert(ss>=cfish->iFishStart && ss<=cfish->iFishEnd);

        		// fill chi
        		// assume width is major axis, else correction:
        		const Real offset = cfish->height[ss] > cfish->width[ss] ? 0.5*M_PI : 0.0;
        		const Real ell_a = (Real)std::max(cfish->height[ss],cfish->width[ss]);
        		//                    const Real dtheta_target = ell_a == 0 ? 2.0*M_PI : 0.25*info.h[0]/ell_a;
        		// maximum distance between two points is 2*ell_a * sin(theta). set this distance to dx/2 -->
        		const Real dtheta_target = ell_a == 0 ? 2.0*M_PI : std::abs(std::asin(0.5*(0.25*info.h_gridpoint)/ell_a));

        		const int Ntheta = (int)std::ceil(2.0*M_PI/dtheta_target) + 1;
        		const Real dtheta = 2.0*M_PI/((Real)Ntheta);

        		for(int tt=0;tt<Ntheta;++tt) {
        			const Real theta = tt*dtheta + offset;
        			// create a surface point
        			// special treatment of tail (width = 0 --> no ellipse, just line)
        			const Real hght = cfish->width[ss] == 0 ? cfish->height[ss]*(2*tt/((Real)Ntheta-1) - 1)
        					: cfish->height[ss]*std::sin(theta);
        			Real myP[3] = {
        					(cfish->rX[ss] + cfish->width[ss]*std::cos(theta)*cfish->norX[ss]),
							(cfish->rY[ss] + cfish->width[ss]*std::cos(theta)*cfish->norY[ss]),
							hght
        			};
        			changeToComputationalFrame(myP);

        			const int iap[3] = {
        					(int)std::floor((myP[0]-org[0])*invh),
							(int)std::floor((myP[1]-org[1])*invh),
							(int)std::floor((myP[2]-org[2])*invh)
        			};
        			// support is two points left, two points right --> Towers Chi will be one point left, one point right, but needs SDF wider
					const int start[3] = {
							std::max(-1, 0 - iap[0] ),
							std::max(-1, 0 - iap[1] ),
							std::max(-1, 0 - iap[2] )
					};
					const int end[3] = {
							std::min(+3, FluidBlock::sizeX - iap[0]),
							std::min(+3, FluidBlock::sizeY - iap[1]),
							std::min(+3, FluidBlock::sizeZ - iap[2])
					};
					const Real myP_distToMidlineSq = cfish->width[ss] == 0 ? std::pow(hght,2) :
							(Real)(std::pow(cfish->width[ss]*std::cos(theta),2) + std::pow(cfish->height[ss]*std::sin(theta),2));

					if(myP_distToMidlineSq<std::numeric_limits<Real>::epsilon()) {
						// if ss==iFishStart or ss==iFishEnd, the surface point lies on the midline --> myP_distToMidlineSq=0.
								// here our in/out criterion fails
						for(int sz=start[2]; sz<end[2];++sz)
						for(int sy=start[1]; sy<end[1];++sy)
						for(int sx=start[0]; sx<end[0];++sx) {
							const int idx[3] = {
									iap[0] + sx,
									iap[1] + sy,
									iap[2] + sz,
							};
							assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
							assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
							assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);

							Real p[3];
							info.pos(p, idx[0],idx[1],idx[2]);
							const Real diff[3] = {p[0]-myP[0], p[1]-myP[1], p[2]-myP[2]};
							const Real distSq = std::pow(diff[0],2) + std::pow(diff[1],2) + std::pow(diff[2],2);
							int closest_s;
							const Real distToMidlineSq = getSmallerDistToMidline(ss, p, closest_s);

							changeFromComputationalFrame(p);
							const Real distPlanar = std::sqrt( std::pow(p[0]-cfish->rX[closest_s],2) +
									std::pow(p[1]-cfish->rY[closest_s],2) );
							const Real distHeight = std::abs(p[2]);
							const Real sign = (distPlanar > cfish->width[closest_s] or distHeight > cfish->height[closest_s]) ? -1.0 : 1.0;

							defblock->chi[idx[2]][idx[1]][idx[0]] =
									(std::abs(defblock->chi[idx[2]][idx[1]][idx[0]]) > distSq) ? sign*distSq
											: defblock->chi[idx[2]][idx[1]][idx[0]];
						}
					} else {
						for(int sz=start[2]; sz<end[2];++sz)
						for(int sy=start[1]; sy<end[1];++sy)
						for(int sx=start[0]; sx<end[0];++sx) {
							const int idx[3] = {
									iap[0] + sx,
									iap[1] + sy,
									iap[2] + sz,
							};
							assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
							assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
							assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);

							Real p[3];
							info.pos(p, idx[0],idx[1],idx[2]);
							const Real diff[3] = {p[0]-myP[0], p[1]-myP[1], p[2]-myP[2]};
							const Real distSq = std::pow(diff[0],2) + std::pow(diff[1],2) + std::pow(diff[2],2);
							int closest_s;
							const Real distToMidlineSq = getSmallerDistToMidline(ss, p, closest_s);
							const Real sign = distToMidlineSq >= myP_distToMidlineSq ? -1.0 : 1.0;
							defblock->chi[idx[2]][idx[1]][idx[0]] =
									(std::abs(defblock->chi[idx[2]][idx[1]][idx[0]]) > distSq) ? sign*distSq
											: defblock->chi[idx[2]][idx[1]][idx[0]];
						}
					}
        		}
        	}
        }
        
        // construct the deformation velocities (P2M with hat function as kernel)
        
        for(int i=0;i<vSegments.size();++i) {
        	for(int ss=vSegments[i].s_range.first;ss<=vSegments[i].s_range.second;++ss) {
        		assert(ss>=0 && ss<=cfish->Nm-1);
        		// P2M udef of a slice at this s
        		const Real myWidth =  (ss < cfish->iFishStart ? cfish->width[ cfish->iFishStart]
									  : (ss > cfish->iFishEnd   ? cfish->width[ cfish->iFishEnd] : cfish->width[ss]));
        		const Real myHeight = (ss < cfish->iFishStart ? cfish->height[cfish->iFishStart]
									  : (ss > cfish->iFishEnd   ? cfish->height[cfish->iFishEnd] : cfish->height[ss]));
        		const Real ds_defGrid = info.h_gridpoint;
        		// towers needs 1dx on each side, smooth needs 2dx --> make it 3 to be nice (and so we can floor!)
        		const Real extension = NPPEXT*info.h_gridpoint; //G tmp changed back to 2

        		const int Nh = std::floor( (myHeight+extension)/ds_defGrid );

        		for(int ih=-Nh;ih<=Nh; ++ih) {
        			const Real offsetH = ih*ds_defGrid;
        			// add an extra extension when width == 0 (to deal with large curvatures near head and/or tail):
        			const Real currentWidth = myWidth== 0 ? extension : myWidth * std::sqrt(1 - std::pow(offsetH/(myHeight+extension),2));
        			const Real actualWidth = (cfish->height[ss] == 0 or std::abs(offsetH)>=cfish->height[ss]) ? 0.0
        										: cfish->width[ss] * std::sqrt(1 - std::pow(offsetH/cfish->height[ss],2));
        			const int Nw = std::floor( (currentWidth+extension)/ds_defGrid); // add extension here to make sure we have it in each direction

        			for(int iw=-Nw;iw<=Nw; ++iw) {
        				const Real offsetW = iw*ds_defGrid;
        				Real xp[3] = {
        						(cfish->rX[ss] + offsetW*cfish->norX[ss]),
								(cfish->rY[ss] + offsetW*cfish->norY[ss]),
								offsetH
        				};
        				changeToComputationalFrame(xp);
        				xp[0] = (xp[0]-org[0])*invh;
        				xp[1] = (xp[1]-org[1])*invh;
        				xp[2] = (xp[2]-org[2])*invh;
        				Real udef[3] = {
        						(cfish->vX[ss] + offsetW*cfish->vNorX[ss]),
								(cfish->vY[ss] + offsetW*cfish->vNorY[ss]),
								0.0
        				};
        				changeVelocityToComputationalFrame(udef);
        				const Real ap[3] = {
        						std::floor((Real)xp[0]),
								std::floor((Real)xp[1]),
								std::floor((Real)xp[2])
        				};
        				const int iap[3] = {
        						(int)ap[0],
								(int)ap[1],
								(int)ap[2]
        				};
        				// now we P2M
						const int start[3] = {
								std::max(0, 0 - iap[0] ),
								std::max(0, 0 - iap[1] ),
								std::max(0, 0 - iap[2] )
        				};
        				const int end[3] = {
        						std::min(+2, FluidBlock::sizeX - iap[0]),
								std::min(+2, FluidBlock::sizeY - iap[1]),
								std::min(+2, FluidBlock::sizeZ - iap[2])
        				};

        				Real wghts[3][2];
        				for(int c=0;c<3;++c) {
        					const Real t[2] = {
        							std::abs((Real)xp[c] - (ap[c]+0)),
									std::abs((Real)xp[c] - (ap[c]+1))
        					};
        					wghts[c][0] = 1.0 - t[0];
        					wghts[c][1] = 1.0 - t[1];
        				}

        				const bool isInside = (std::abs(offsetW) < actualWidth) && (std::abs(offsetH) < cfish->height[ss]);
        				for(int sz=start[2]; sz<end[2];++sz) {
        					const Real wz = wghts[2][sz];
        					for(int sy=start[1];sy<end[1];++sy) {
        						const Real wywz = wz*wghts[1][sy];
        						for(int sx=start[0];sx<end[0];++sx) {
        							const Real wxwywz = wywz*wghts[0][sx];
        							assert(wxwywz>=0 && wxwywz<=1);
        							const int idx[3] = {
        									iap[0] + sx,
											iap[1] + sy,
											iap[2] + sz,
        							};
        							assert(idx[0]>=0 && idx[0]<FluidBlock::sizeX);
        							assert(idx[1]>=0 && idx[1]<FluidBlock::sizeY);
        							assert(idx[2]>=0 && idx[2]<FluidBlock::sizeZ);
        							defblock->udef[idx[2]][idx[1]][idx[0]][0] += wxwywz*udef[0];
        							defblock->udef[idx[2]][idx[1]][idx[0]][1] += wxwywz*udef[1];
        							defblock->udef[idx[2]][idx[1]][idx[0]][2] += wxwywz*udef[2];
        							b(idx[0],idx[1],idx[2]).tmpU += wxwywz;
        							// set sign for all interior points
									if( (std::abs(defblock->chi[idx[2]][idx[1]][idx[0]] + 1) < 5*std::numeric_limits<Real>::epsilon()) && isInside)
										defblock->chi[idx[2]][idx[1]][idx[0]] = 1.0;
        						}
        					}
        				}
        			}
        		}
        	}
        }
        
        // finalize signed distance function in tmpU
        {
            for(int iz=0; iz<FluidBlock::sizeZ; iz++)
			for(int iy=0; iy<FluidBlock::sizeY; iy++)
			for(int ix=0; ix<FluidBlock::sizeX; ix++) {
				const Real normfac = b(ix,iy,iz).tmpU > numeric_limits<Real>::epsilon() ? b(ix,iy,iz).tmpU : 1.;
				defblock->udef[iz][iy][ix][0] /= normfac;
				defblock->udef[iz][iy][ix][1] /= normfac;
				defblock->udef[iz][iy][ix][2] /= normfac;
				// change from signed squared distance function to normal sdf
				b(ix,iy,iz).tmpU = defblock->chi[iz][iy][ix] > 0 ?
						 	 sqrt( defblock->chi[iz][iy][ix]) :
							-sqrt(-defblock->chi[iz][iy][ix]);
				b(ix,iy,iz).tmpV = defblock->udef[iz][iy][ix][0];
				b(ix,iy,iz).tmpW = defblock->udef[iz][iy][ix][1];
			}
        }
	}
};

struct PutFishOnBlocks_Finalize : public GenericLabOperator
{
    Real t;
	int stencil_start[3], stencil_end[3];
	array<Real,4>* const momenta;
    surfaceBlocks* const surface;
    std::map<int,ObstacleBlock*>* const obstacleBlocks;


	//PutFishOnBlocks_Finalize finalize(obstacleBlocks,dataPerThread[tid],tmp,blockID);

    PutFishOnBlocks_Finalize(map<int,ObstacleBlock*>* const obstacleBlocks, //to write chi
    						 surfaceBlocks* const surface, 					//to write gradChi
							 array<Real,4>* const momenta)     			//compute CM
	: t(0), momenta(momenta),surface(surface), obstacleBlocks(obstacleBlocks)
	{
    	stencil = StencilInfo(-1,-1,-1, 2,2,2, false, 1, 5);
		stencil_start[0] = stencil_start[1] = stencil_start[2] = -1;
		stencil_end[0] = stencil_end[1] = stencil_end[2] = +2;
	}
#if 1
    inline Real sign(const Real& val) const {
    	return (0. < val) - (val < 0.);
    }
   template <typename Lab, typename BlockType>
	void operator()(const Lab& lab, const BlockInfo& info, BlockType& b) const
	{
		if(obstacleBlocks->find(info.blockID) == obstacleBlocks->end()) return;
		ObstacleBlock* const defblock = obstacleBlocks->find(info.blockID)->second;

		const Real eps = std::numeric_limits<Real>::epsilon();
		const Real h = info.h_gridpoint;
		const Real fac1 = .5/h;
		const Real fac2 = .5/(std::sqrt(2.)*h);
		const Real fac[9] = {fac1,fac1,fac1,fac2,fac2,fac2,fac2,fac2,fac2};

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++) {
			Real p[3];
			info.pos(p, ix,iy,iz);
			if (lab(ix,iy,iz).tmpU > +2*h || lab(ix,iy,iz).tmpU < -2*h) {
				const Real H = lab(ix,iy,iz).tmpU > 0 ? 1.0 : 0.0;
				(*momenta)[0] += H;
				(*momenta)[1] += p[0]*H;
				(*momenta)[2] += p[1]*H;
				(*momenta)[3] += p[2]*H;
				b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
				defblock->chi[iz][iy][ix] = H;
				continue;
			}
			const Real Up[9] = {
				lab(ix+1,iy  ,iz  ).tmpU, lab(ix  ,iy+1,iz  ).tmpU, lab(ix  ,iy  ,iz+1).tmpU,
				lab(ix  ,iy+1,iz+1).tmpU, lab(ix+1,iy  ,iz+1).tmpU, lab(ix+1,iy+1,iz  ).tmpU,
				lab(ix  ,iy+1,iz-1).tmpU, lab(ix+1,iy  ,iz-1).tmpU, lab(ix+1,iy-1,iz  ).tmpU
			};
			const Real Um[9] = {
				lab(ix-1,iy  ,iz  ).tmpU, lab(ix  ,iy-1,iz  ).tmpU, lab(ix  ,iy  ,iz-1).tmpU,
				lab(ix  ,iy-1,iz-1).tmpU, lab(ix-1,iy  ,iz-1).tmpU, lab(ix-1,iy-1,iz  ).tmpU,
				lab(ix  ,iy-1,iz+1).tmpU, lab(ix-1,iy  ,iz+1).tmpU, lab(ix-1,iy+1,iz  ).tmpU
			};
			Real Ip[9],Im[9];
			for (int i=0; i<9; i++) {
				Ip[i] = std::max(Up[i],(Real)0.);
				Im[i] = std::max(Um[i],(Real)0.);
			}
			Real gradU[9], gradI[9];
			for (int i=0; i<9; i++) {
				gradU[i] = fac[i]*(Up[i]-Um[i]);
				gradI[i] = fac[i]*(Ip[i]-Im[i]);
			}
			Real Hp[3], Hm[3];
			for(int i=0; i<3; i++) {
				Hp[i] = (Up[i]> h) ? h : (
						(Up[i]<-h) ? 0 :
						.5*h+(Up[i]-fac1*sign(Up[i])*Up[i]*Up[i]));
				Hm[i] = (Um[i]> h) ? h : (
						(Um[i]<-h) ? 0 :
						.5*h+(Um[i]-fac1*sign(Um[i])*Um[i]*Um[i]));
			}
			const Real gradH[3] = {.5*(Hp[0]-Hm[0]), .5*(Hp[1]-Hm[1]), .5*(Hp[2]-Hm[2])};
			Real gradUU[3], gradUI[3], gradUH[3];
			for (int i=0; i<3; i++) {
				gradUU[i] = gradU[i]*gradU[i] + gradU[i+3]*gradU[i+3] + gradU[i+6]*gradU[i+6];
				gradUI[i] = gradU[i]*gradI[i] + gradU[i+3]*gradI[i+3] + gradU[i+6]*gradI[i+6];
				gradUH[i] = gradU[i]*gradH[i];
			}
			for (int i=0; i<3; i++)  gradUU[i] = max(gradUU[i], eps);
			const Real FDD = h*(gradUH[0] + gradUH[1] + gradUH[2])/gradUU[0];
			const Real FDH = 1/3. * (gradUI[0]/gradUU[0]+gradUI[1]/gradUU[1]+gradUI[2]/gradUU[2]);


			if (FDD>eps) {
				const Real dchidx = -FDD*gradU[0];
				const Real dchidy = -FDD*gradU[1];
				const Real dchidz = -FDD*gradU[2];
				surface->add(info.blockID, ix, iy, iz, dchidx, dchidy, dchidz, FDD);
			}
#ifndef NDEBUG
			if(FDH<0 || FDH>1) printf("invalid H?: %9.9e %9.9e %9.9e: %9.9e\n",x,y,z,FDH);
#endif
			(*momenta)[0] += FDH;
			(*momenta)[1] += p[0]*FDH;
			(*momenta)[2] += p[1]*FDH;
			(*momenta)[3] += p[2]*FDH;
			defblock->chi[iz][iy][ix] = FDH;
			b(ix,iy,iz).chi = std::max(FDH, b(ix,iy,iz).chi);
		}
	}
#else
    template <typename Lab, typename BlockType>
	void operator()(const Lab& lab, const BlockInfo& info, BlockType& b)
	{
		if(obstacleBlocks->find(info.blockID) == obstacleBlocks->end()) return;
		ObstacleBlock* const defblock = obstacleBlocks->find(info.blockID)->second;
		const Real eps = std::numeric_limits<Real>::epsilon();
		const Real h = info.h_gridpoint;
		const Real invh2 = 1./(h*h);
		const Real inv2h = .5/h;

		for(int iz=0; iz<FluidBlock::sizeZ; iz++)
		for(int iy=0; iy<FluidBlock::sizeY; iy++)
		for(int ix=0; ix<FluidBlock::sizeX; ix++) {
			Real p[3];
			info.pos(p, ix,iy,iz);
			if (lab(ix,iy,iz).tmpU >= +2*h || lab(ix,iy,iz).tmpU <= -2*h) {
				const Real H = lab(ix,iy,iz).tmpU > 0 ? 1.0 : 0.0;
				(*momenta)[0] += H;
				(*momenta)[1] += p[0]*H;
				(*momenta)[2] += p[1]*H;
				(*momenta)[3] += p[2]*H;
				//printf("The fuck %f %f %f %f\n",H,p[0],p[1],p[2]);
				//printf("The fuck %f %f %f %f\n",normfac,b(ix,iy,iz).tmpU,b(ix,iy,iz).tmpV,b(ix,iy,iz).tmpW);
				b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
				defblock->chi[iz][iy][ix] = H;
				continue;
			}

			const Real distPx = lab(ix+1,iy,iz).tmpU;
			const Real distMx = lab(ix-1,iy,iz).tmpU;
			const Real distPy = lab(ix,iy+1,iz).tmpU;
			const Real distMy = lab(ix,iy-1,iz).tmpU;
			const Real distPz = lab(ix,iy,iz+1).tmpU;
			const Real distMz = lab(ix,iy,iz-1).tmpU;

			const Real IplusX = distPx < 0 ? 0 : distPx;
			const Real IminuX = distMx < 0 ? 0 : distMx;
			const Real IplusY = distPy < 0 ? 0 : distPy;
			const Real IminuY = distMy < 0 ? 0 : distMy;
			const Real IplusZ = distPz < 0 ? 0 : distPz;
			const Real IminuZ = distMz < 0 ? 0 : distMz;

			const Real HplusX = distPx == 0 ? 0.5 : (distPx < 0 ? 0 : 1);
			const Real HminuX = distMx == 0 ? 0.5 : (distMx < 0 ? 0 : 1);
			const Real HplusY = distPy == 0 ? 0.5 : (distPy < 0 ? 0 : 1);
			const Real HminuY = distMy == 0 ? 0.5 : (distMy < 0 ? 0 : 1);
			const Real HplusZ = distPz == 0 ? 0.5 : (distPz < 0 ? 0 : 1);
			const Real HminuZ = distMz == 0 ? 0.5 : (distMz < 0 ? 0 : 1);

			// gradU
			const Real gradUX = inv2h * (distPx - distMx);
			const Real gradUY = inv2h * (distPy - distMy);
			const Real gradUZ = inv2h * (distPz - distMz);
			const Real gradUSq = gradUX*gradUX + gradUY*gradUY + gradUZ*gradUZ;

			// gradI: first primitive of H(x): I(x) = int_0^x H(y) dy
			const Real gradIX = inv2h * (IplusX - IminuX);
			const Real gradIY = inv2h * (IplusY - IminuY);
			const Real gradIZ = inv2h * (IplusZ - IminuZ);
			const Real numH = gradIX*gradUX + gradIY*gradUY + gradIZ*gradUZ;

			const Real gradHX = inv2h * (HplusX - HminuX);
			const Real gradHY = inv2h * (HplusY - HminuY);
			const Real gradHZ = inv2h * (HplusZ - HminuZ);
			const Real numD = gradHX*gradUX + gradHY*gradUY + gradHZ*gradUZ;

			const Real Delta = std::abs(gradUSq) < eps ? numD : numD/gradUSq;
			const Real H     = std::abs(gradUSq) < eps ? numH : numH/gradUSq;

			if (Delta>1e-6) {
				const Real dchidx = -Delta*gradUX;
				const Real dchidy = -Delta*gradUY;
				const Real dchidz = -Delta*gradUZ;
				surface->add(info.blockID, ix, iy, iz, dchidx, dchidy, dchidz, Delta);
			}
			(*momenta)[0] += H;
			(*momenta)[1] += p[0]*H;
			(*momenta)[2] += p[1]*H;
			(*momenta)[3] += p[2]*H;
			defblock->chi[iz][iy][ix] = H;
			b(ix,iy,iz).chi = std::max(H, b(ix,iy,iz).chi);
		}
	}
#endif
};

void IF3D_CarlingFishOperator::create(const int step_id,const Real time, const Real dt, const Real *Uinf)
{
    // STRATEGY
    // we need some things already
    // - the internal angle at the previous timestep, obtained from integrating the actual def velocities
	// 						 (not the imposed deformation velocies, because they dont have zero ang mom)
    // - the internal angular velocity at previous timestep
    
    // 1. create midline
    // 2. integrate to find CoM, angular velocity, etc
    // 3. shift midline to CoM frame: zero internal linear momentum and angular momentum
    
    // 4. split the fish into segments (according to s)
    // 5. rotate the segments to computational frame (comp CoM and angle)
    // 6. for each Block in the domain, find those segments that intersect it
    // 7. for each of those blocks, allocate an ObstacleBlock
    
    // 8. put the 3D shape on the grid: SDF-P2M for sdf, normal P2M for udef
    // 9. create the Chi out of the SDF. In same sweep, compute the actual CoM
    // 10. compute all shit: linear momentum, angular momentum etc.
    // 11. correct deformation velocity to nullify momenta for the final discrete representation
    
    const int Nsegments = NPPSEG;
    const int Nextension = 4*NPPEXT;// up to 3dx on each side (to get proper interpolation up to 2dx)
    const Real dx_extension = 0.25*vInfo[0].h_gridpoint;
    const Real target_ds = vInfo[0].h_gridpoint/TGTPPB;
    const Real target_Nm = length/target_ds;
    const int Nm = NPPSEG*(int)std::ceil(target_Nm/NPPSEG) + 1;
    assert((Nm-1)%Nsegments==0);
    if (bCorrectTrajectory) {
        Real velx_tot = Uinf[0] - transVel[0];
        Real vely_tot = Uinf[1] - transVel[1];
        Real AngDiff  = atan2(vely_tot,velx_tot);
        adjTh = (1.-dt) * adjTh + dt * AngDiff;
        const Real B = (AngDiff*angVel[2]>0) ? 0.25/M_PI : 0;
        const Real PID = .5*adjTh +B*AngDiff*fabs(angVel[2]);
        myFish->_correctTrajectory(PID, time, dt);
    }
    // 1.
    myFish->computeMidline(time);
    
    // 2. & 3.
    volume_internal = myFish->integrateLinearMomentum(CoM_internal, vCoM_internal);
    assert(volume_internal > std::numeric_limits<Real>::epsilon());
    myFish->changeToCoMFrameLinear(CoM_internal, vCoM_internal);
    
    angvel_internal_prev = angvel_internal;
    myFish->integrateAngularMomentum(angvel_internal);
    J_internal = myFish->J;
    // update theta now with new angvel info
    //theta_internal -= 0.5*sim_dt*(angvel_internal+angvel_internal_prev);//negative: we subtracted this angvel
    myFish->changeToCoMFrameAngular(theta_internal, angvel_internal);

#ifndef NDEBUG
    {
        Real dummy_CoM_internal[2], dummy_vCoM_internal[2], dummy_angvel_internal;
        // check that things are zero
        const Real volume_internal_check = myFish->integrateLinearMomentum(dummy_CoM_internal,dummy_vCoM_internal);
        myFish->integrateAngularMomentum(dummy_angvel_internal);
        
        assert(std::abs(dummy_CoM_internal[0])<10*std::numeric_limits<Real>::epsilon());
        assert(std::abs(dummy_CoM_internal[1])<10*std::numeric_limits<Real>::epsilon());
        assert(std::abs(myFish->linMom[0])<10*std::numeric_limits<Real>::epsilon());
        assert(std::abs(myFish->linMom[1])<10*std::numeric_limits<Real>::epsilon());
        assert(std::abs(myFish->angMom)<10*std::numeric_limits<Real>::epsilon());
        assert(std::abs(volume_internal - volume_internal_check) < 10*std::numeric_limits<Real>::epsilon());
    }
#endif
    
    // 4.
    std::vector<VolumeSegment_OBB> vSegments(Nsegments);
#pragma omp parallel for
	for(int i=0;i<Nsegments;++i) {
		const int next_idx = (i+1)*(Nm-1)/Nsegments;
		const int idx = i * (Nm-1)/Nsegments;
		// find bounding box based on this
		Real bbox[3][2] = {{1e9, -1e9}, {1e9, -1e9}, {1e9, -1e9}};
		for(int ss=idx; ss<=next_idx; ++ss) {
			const Real xBnd[2] = {myFish->rX[ss] - myFish->norX[ss]*myFish->width[ss],
									myFish->rX[ss] + myFish->norX[ss]*myFish->width[ss]};
			const Real yBnd[2] = {myFish->rY[ss] - myFish->norY[ss]*myFish->width[ss],
									myFish->rY[ss] + myFish->norY[ss]*myFish->width[ss]};
			const Real zBnd[2] = {-myFish->height[ss], +myFish->height[ss]};
			bbox[0][0] = std::min({bbox[0][0],xBnd[0],xBnd[1]});
			bbox[0][1] = std::max({bbox[0][1],xBnd[0],xBnd[1]});
			bbox[1][0] = std::min({bbox[1][0],yBnd[0],yBnd[1]});
			bbox[1][1] = std::max({bbox[1][1],yBnd[0],yBnd[1]});
			bbox[2][0] = std::min({bbox[2][0],zBnd[0],zBnd[1]});
			bbox[2][1] = std::max({bbox[2][1],zBnd[0],zBnd[1]});
		}
		vSegments[i].prepare(std::make_pair(idx, next_idx), bbox); //create a new segment
		vSegments[i].changeToComputationalFrame(position,quaternion);
	}

    /*
    std::vector<VolumeSegment_OBB> vSegments;
    vSegments.reserve(Nsegments);
    {
        for(int i=0;i<Nsegments;++i) {
            const int next_idx = (i+1)*(Nm-1)/Nsegments;
            const int idx = i * (Nm-1)/Nsegments;

            for(int ss=idx; ss<=next_idx; ++ss) {
                const Real xBnd[2] = {myFish->rX[ss] - myFish->norX[ss]*myFish->width[ss],
                						myFish->rX[ss] + myFish->norX[ss]*myFish->width[ss]};
                const Real yBnd[2] = {myFish->rY[ss] - myFish->norY[ss]*myFish->width[ss],
                						myFish->rY[ss] + myFish->norY[ss]*myFish->width[ss]};
                const Real zBnd[2] = {-myFish->height[ss], +myFish->height[ss]};
                bbox[0][0] = std::min({bbox[0][0],xBnd[0],xBnd[1]});
                bbox[0][1] = std::max({bbox[0][1],xBnd[0],xBnd[1]});
                bbox[1][0] = std::min({bbox[1][0],yBnd[0],yBnd[1]});
                bbox[1][1] = std::max({bbox[1][1],yBnd[0],yBnd[1]});
                bbox[2][0] = std::min({bbox[2][0],zBnd[0],zBnd[1]});
                bbox[2][1] = std::max({bbox[2][1],zBnd[0],zBnd[1]});
            }
            
            VolumeSegment_OBB volumeSegment(std::make_pair(idx, next_idx), bbox); //create a new segment
            vSegments.push_back(volumeSegment);
        }
    }
    assert(vSegments.size()==Nsegments);

    // 5.
#pragma omp parallel for
    for(int i=0;i<Nsegments;++i)
        vSegments[i].changeToComputationalFrame(position,quaternion);
    */

    // clear deformation velocities
    for(auto & entry : obstacleBlocks)
        delete entry.second;
    obstacleBlocks.clear();

    // 6. & 7.    
    std::map<int, std::vector<VolumeSegment_OBB>> segmentsPerBlock;
    {
        for(int i=0;i<vInfo.size();++i) {
            const BlockInfo & info = vInfo[i];
            Real pStart[3], pEnd[3];
            info.pos(pStart, 0, 0, 0);
            info.pos(pEnd, FluidBlock::sizeX-1, FluidBlock::sizeY-1, FluidBlock::sizeZ-1);
            const Real safe_distance = 2.0*info.h_gridpoint; // two points on each side

            for(int s=0;s<Nsegments;++s)
                if(vSegments[s].isIntersectingWithAABB(pStart,pEnd,safe_distance))
                    segmentsPerBlock[info.blockID].push_back(vSegments[s]);
            
            // allocate new blocks if necessary
            if(segmentsPerBlock.find(info.blockID) != segmentsPerBlock.end()) {
                assert(obstacleBlocks.find(info.blockID) == obstacleBlocks.end());
                obstacleBlocks[info.blockID] = new ObstacleBlock;
                obstacleBlocks[info.blockID]->clear();
            }
        }
    }
    
    //assert(not segmentsPerBlock.empty()); //killed this assert: distributed fish
    assert(segmentsPerBlock.size() == obstacleBlocks.size());

    // 8.
    {
#pragma omp parallel
    	{
        	PutFishOnBlocks putfish(myFish, position, quaternion);

#pragma omp for schedule(static)
    		for(int i=0; i<vInfo.size(); i++) {
    			BlockInfo info = vInfo[i];
    			auto pos = segmentsPerBlock.find(info.blockID);
    			FluidBlock& b = *(FluidBlock*)info.ptrBlock;

    			//tmpU will contain SDF: neg outside positive inside
    			if(pos == segmentsPerBlock.end()) {
    				for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
					for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
						b(ix,iy,iz).tmpU = -1.; //-1 here to avoid gremlins at blocks' boundaries
    			} else {
    				for(int iz=0; iz<FluidBlock::sizeZ; ++iz)
					for(int iy=0; iy<FluidBlock::sizeY; ++iy)
					for(int ix=0; ix<FluidBlock::sizeX; ++ix)
						b(ix,iy,iz).tmpU = 0.; //this will be accessed with plus equal

					assert(obstacleBlocks.find(info.blockID) != obstacleBlocks.end());
					ObstacleBlock* const defblock = obstacleBlocks.find(info.blockID)->second;
					putfish(info, b, defblock, pos->second);
    			}
    		}
    	}
    }

    // 9. & 10. & 11.
    {
    	const int nthreads = omp_get_max_threads();
    	vector<surfaceBlocks> dataPerThread(nthreads);
    	vector<array<Real,4>> momenta(nthreads);
    	vector<PutFishOnBlocks_Finalize*> finalize;
    	for(int i = 0; i < nthreads; ++i) {
    		PutFishOnBlocks_Finalize* tmp = new PutFishOnBlocks_Finalize(&obstacleBlocks,&dataPerThread[i],&momenta[i]);
    		finalize.push_back(tmp);
    	}
    	compute(finalize);

    	double sumX[4] = {0,0,0,0};
    	double totX[4] = {0,0,0,0};
    	for(int i=0; i<nthreads; i++) {
    		sumX[0] += momenta[i][0];
    		sumX[1] += momenta[i][1];
    		sumX[2] += momenta[i][2];
    		sumX[3] += momenta[i][3];
    	}

    	MPI::COMM_WORLD.Allreduce(sumX, totX, 4, MPI::DOUBLE, MPI::SUM);

    	surfData.finalizeOnGrid(dataPerThread);

        assert(totX[0]>std::numeric_limits<double>::epsilon());
        CoM_interpolated[0]=totX[1]/totX[0];
        CoM_interpolated[1]=totX[2]/totX[0];
        CoM_interpolated[2]=totX[3]/totX[0];

        _makeDefVelocitiesMomentumFree(CoM_interpolated);
    }
}

void IF3D_CarlingFishOperator::update(const int stepID, const Real t, const Real dt, const Real *Uinf)
{
    // synchronize internal time
    sim_time = t + dt;
    sim_dt = dt;
    // update position and angles
    IF3D_ObstacleOperator::update(stepID,t, dt, Uinf);
    // negative: we subtracted this angvel
    theta_internal -= sim_dt*angvel_internal;
}

void IF3D_CarlingFishOperator::_parseArguments(ArgumentParser & parser)
{
	IF3D_ObstacleOperator::_parseArguments(parser);
	parser.set_strict_mode();
	Tperiod = parser("-T").asDouble();
    parser.unset_strict_mode();
    nActions = parser("-nActions").asInt(0);
    GoalDX = parser("-GoalDX").asDouble(0.0);
	phaseShift = parser("-phi").asDouble(0.0);
    Tstartlearn = parser("-Tstartlearn").asDouble(1e6);
	bCorrectTrajectory = parser("-Correct").asBool(false);
	/*
    randomStart = parser("-randomStart").asBool(false);
    if (randomStart) {
    	printf("Random start\n");
    	std::random_device rd;
    	std::mt19937 gen(rd());
    	std::uniform_real_distribution<Real> dis(-1.,1.);
    	const Real rseed1 = .8*length*dis(gen);
    	const Real rseed2 = .2*length*dis(gen);
    	const Real rseed3 = .1* M_PI *dis(gen);
    	position[0] += rseed1/sqrt(2.)-rseed2/sqrt(2.);
    	position[1] += rseed1/sqrt(2.)+rseed2/sqrt(2.);
    	ext_pos[0] = position[0];
    	ext_pos[1] = position[1];
    	angle += rseed3;
    }
    sr->updateInstant(position[0], ext_pos[0], position[1], ext_pos[1], angle, 0., 0., 0.);
    //TODO state and reward:
    sr->t_next_comm = Tstartlearn - 1/2.; //i want to reset time-averages before first actual comm
    bool bForgiving = parser("-easyFailBox").asBool(false);
    sr->bForgiving = bForgiving;
    sr->GoalDX = GoalDX;
    sr->thExp = angle;

    randomActions = parser("-randomActions").asBool(false);
    if (randomActions) printf("Fish doing random turns\n");
    useLoadedActions = parser("-useLoadedActions").asBool(false);
    if (useLoadedActions) {
        Real dummy_time;
        vector<Real> action(nActions);
        ifstream in("orders_1.txt"); //FUCKING TODO NEED TO USE SOME POINTERS IN THIS SHIT
        std::string line;
        if(in.good()) {
            while (getline(in, line)) {
                istringstream line_in(line);
                line_in >> dummy_time;
                line_in >> action[0];
                if(nActions==2) line_in >> action[1];
                //i want to do pop back later:
                loadedActions.insert(loadedActions.begin(),action);
            }
        } else { printf("Could not load actions from file orders_1.txt\n"); abort(); }
        in.close();
    }
    */
}

/*
void IF3D_StefanLearnTurnOperator::execute(Communicator * comm, const int iAgent, const Real time)
{
    if (time < Tstartlearn) {
        sr->resetAverage();
        sr->t_next_comm = Tstartlearn;

        //TMP: first rnd action is going to be taken after a while
        if (randomActions) sr->t_next_comm = Tstartlearn+6.;

        return;
    }

    if (not bInteractive) {
        if (not randomActions) { sr->t_next_comm=1e6; return; }
        //we might decide to pause turning, then just pick a pause counter
        if (nPauseActions-- > 0) {
            vector<Real> raction(1,0.);
            myFish->execute(time, sr->t_next_comm, raction);
            sr->t_next_comm += .5*myFish->Tperiod;
            printf("pausing from turning at time %f, will still pause %d turns.
            		Next turn at time %f\n",time,nPauseActions,sr->t_next_comm);
            return;
        }
        vector<Real> raction(1,signLastTurn);
        myFish->execute(time, sr->t_next_comm, raction);
        sr->t_next_comm += .5*myFish->Tperiod;
        printf("turning at time %f with modifier %f, turn counter is %d.
        		Next turn at time %f\n",time,signLastTurn,nTurnActions,sr->t_next_comm);

        if (++nTurnActions >= 4) {
            nPauseActions = 10.;
            nTurnActions = 0;
            signLastTurn *= -1.;
        }

    } else if (useLoadedActions) {

        vector<Real> actions(nActions);
        if (loadedActions.size()>1) {
            actions = loadedActions.back();
            loadedActions.pop_back();
        } //else zero actions
        myFish->execute(time, sr->t_next_comm, actions);

        old_curv = new_curv;
        new_curv = actions[0];
        if(nActions==2) {
            new_Tp = actions[1];
            sr->t_next_comm += .5*myFish->l_Tp;
        } else if (nActions==1) {
            sr->t_next_comm += .5*myFish->Tperiod;
        }
        sr->resetAverage();

    } else {

        const Real relT= fmod(time,1.); //1 is Tperiod
#ifdef _NOVISION_
        const int nStates = (nActions==1) ? 20+ 8*NpLatLine : 25+  8*NpLatLine;
#else
        const int nStates = (nActions==1) ? 20+10*NpLatLine : 25+ 10*NpLatLine;
#endif
        vector<Real> state(nStates), actions(nActions);

        int k(0);
        state[k++] = sr->Xrel - GoalDX;
        state[k++] = sr->Yrel;
        state[k++] = sr->RelAng;
        state[k++] = relT;
        state[k++] = new_curv;
        state[k++] = old_curv;

        if(nActions==2) { //this is for backwards compatibility
            state[k++] = new_Tp;
                        //2.*M_PI*((time-time0)/l_Tp +timeshift -rS[i]/length) + M_PI*phaseShift
            Real Fshift = 2.*((-myFish->time0)/myFish->l_Tp +myFish->timeshift)+myFish->phaseShift;
            Fshift = fmod(Fshift,2.0);
            state[k++] = (Fshift<0) ? 2.+Fshift : Fshift;
            state[k++] = sr->VX;
            state[k++] = sr->VY;
            state[k++] = sr->AV;
        }

        state[k++] = sr->Dist;
        state[k++] = sr->Quad;
        state[k++] = sr->VxAvg;
        state[k++] = sr->VyAvg;
        state[k++] = sr->AvAvg;
        state[k++] = sr->Pout;
        state[k++] = sr->defPower;
        state[k++] = sr->EffPDef;
        state[k++] = sr->PoutBnd;
        state[k++] = sr->defPowerBnd;
        state[k++] = sr->EffPDefBnd;
        state[k++] = sr->Pthrust;
        state[k++] = sr->Pdrag;
        state[k++] = sr->ToD;

        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelNAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelTAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelNBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->VelTBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FPAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FVAbove[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FPBelow[j];
        for (int j=0; j<NpLatLine; j++) state[k++] = sr->FVBelow[j];
        #ifndef _NOVISION_
        for (int j=0; j<2*NpLatLine; j++) state[k++] = sr->raySight[j];
        #endif
        const Real reward = (sr->info==2) ? -10 : sr->EffPDefBnd;
        comm->sendState(iAgent-1, sr->info, state, reward); //TODO
        fflush(0);
        if (sr->info==2) return;

        sr->info = 0;

        comm->recvAction(actions);
        myFish->execute(time, sr->t_next_comm, actions);

        old_curv = new_curv;
        new_curv = actions[0];
        if(nActions==2) {
            new_Tp = actions[1];
            sr->t_next_comm += .5*myFish->l_Tp;
        } else if (nActions==1) {
            sr->t_next_comm += .5*myFish->Tperiod;
        }

        #ifndef TRAINING
        ofstream filedrag;
        filedrag.open(("orders_"+to_string(iAgent)+".txt").c_str(), ios::app);
        filedrag<<time<<" "<<new_curv;
        if(nActions==2)
            filedrag<<" "<<new_Tp;
        filedrag<<endl;
        filedrag.close();
        #endif //TRAINING

        sr->resetAverage();
    }
}
*/

void IF3D_CarlingFishOperator::getCenterOfMass(Real CM[3]) const
{
	// return computation CoM, not the one were advecting
	CM[0]=CoM_interpolated[0];
	CM[1]=CoM_interpolated[1];
	CM[2]=CoM_interpolated[2];
}

void IF3D_CarlingFishOperator::save(const int stepID, const Real t, string filename)
{
    //assert(std::abs(t-sim_time)<std::numeric_limits<Real>::epsilon());
    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<Real>::digits10 + 1);
    savestream.open(filename + ".txt");

    savestream<<t<<"\t"<<sim_dt<<std::endl;
    savestream<<position[0]<<"\t"<<position[1]<<"\t"<<position[2]<<std::endl;
    savestream<<quaternion[0]<<"\t"<<quaternion[1]<<"\t"<<quaternion[2]<<"\t"<<quaternion[3]<<std::endl;
    savestream<<transVel[0]<<"\t"<<transVel[1]<<"\t"<<transVel[2]<<std::endl;
    savestream<<angVel[0]<<"\t"<<angVel[1]<<"\t"<<angVel[2]<<std::endl;
    savestream<<theta_internal<<"\t"<<angvel_internal<<"\t"<<adjTh<<std::endl;
    savestream.close();
    
}

void IF3D_CarlingFishOperator::restart(const Real t, string filename)
{
    std::ifstream restartstream;
    restartstream.open(filename+".txt");
    restartstream >> sim_time >> sim_dt;
    assert(std::abs(sim_time-t) < std::numeric_limits<Real>::epsilon());
    restartstream >> position[0] >> position[1] >> position[2];
    restartstream >> quaternion[0] >> quaternion[1] >> quaternion[2] >> quaternion[3];
    restartstream >> transVel[0] >> transVel[1] >> transVel[2];
    restartstream >> angVel[0] >> angVel[1] >> angVel[2];
    restartstream >> theta_internal >> angvel_internal >> adjTh;
    restartstream.close();
    
	std::cout<<"RESTARTED FISH: "<<std::endl;
	std::cout<<"TIME, DT: "<<sim_time<<" "<<sim_dt<<std::endl;
	std::cout<<"POS: "<<position[0]<<" "<<position[1]<<" "<<position[2]<<std::endl;
	std::cout<<"ANGLE: "<<quaternion[0]<<" "<<quaternion[1]<<" "<<quaternion[2]<<" "<<quaternion[3]<<std::endl;
	std::cout<<"TVEL: "<<transVel[0]<<" "<<transVel[1]<<" "<<transVel[2]<<std::endl;
	std::cout<<"AVEL: "<<angVel[0]<<" "<<angVel[1]<<" "<<angVel[2]<<std::endl;
	std::cout<<"INTERN: "<<theta_internal<<" "<<angvel_internal<<std::endl;
}

/*
void IF3D_CarlingFishOperator::writeToFile(const int step_id, const Real t, std::string filename)
{
    std::string fname = (filename==std::string()) ? "fish" : filename;

    std::ofstream savestream;
    savestream.setf(std::ios::scientific);
    savestream.precision(std::numeric_limits<Real>::digits10 + 1);
    
    savestream.open(fname+"_interpolated.dat", ios::app | ios::out);
    if(step_id==0)
    {
        savestream << "step\t";
        savestream << "time\t";
        savestream << "volume\t";        
        savestream << "CoM[0]\t";
        savestream << "CoM[1]\t";
        savestream << "CoM[2]\t";
        savestream << "linMom[0]\t";
        savestream << "linMom[1]\t";
        savestream << "linMom[2]\t";
        savestream << "angMom[0]\t";
        savestream << "angMom[1]\t";
        savestream << "angMom[2]" << std::endl;
    }    
        
    savestream << step_id << "\t";
    savestream << sim_time << "\t";
    savestream << object_ongrid.volume << "\t";
    savestream << CoM_interpolated[0] << "\t";
    savestream << CoM_interpolated[1] << "\t";
    savestream << CoM_interpolated[2] << "\t";
    savestream << object_ongrid.linearMomentum[0] << "\t";
    savestream << object_ongrid.linearMomentum[1] << "\t";
    savestream << object_ongrid.linearMomentum[2] << "\t";
    savestream << object_ongrid.angularMomentum[0] << "\t";
    savestream << object_ongrid.angularMomentum[1] << "\t";
    savestream << object_ongrid.angularMomentum[2] << "\t";
    savestream << object_ongrid.J[2] << std::endl;
    savestream.close();
    
    savestream.open(fname+"_internal.dat", ios::app | ios::out);
    if(step_id==0)
    {
        savestream << "step\t";
        savestream << "time\t";
        savestream << "volume\t";
        savestream << "CoM[0]\t";
        savestream << "CoM[1]\t";
        savestream << "linMom[0]\t";
        savestream << "linMom[1]\t";
        savestream << "angMom\t";
        savestream << "theta\t";
        savestream << "angvel" << std::endl;
    }
    
    savestream << step_id << "\t";
    savestream << sim_time << "\t";
    savestream << volume_internal << "\t";
    savestream << CoM_internal[0] << "\t";
    savestream << CoM_internal[1] << "\t";
    savestream << vCoM_internal[0]*volume_internal << "\t";
    savestream << vCoM_internal[1]*volume_internal << "\t";
    savestream << angvel_internal*J_internal << "\t";
    savestream << theta_internal << "\t";
    savestream << angvel_internal << std::endl;
    savestream.close();
    
    savestream.open(fname+"_computation.dat", ios::app | ios::out);
    if(step_id==0)
    {
        savestream << "step\t";
        savestream << "time\t";
        savestream << "pos[0]\t";
        savestream << "pos[1]\t";
        savestream << "pos[2]\t";
        savestream << "quat[0]\t";
        savestream << "quat[1]\t";
        savestream << "quat[2]\t";
        savestream << "quat[3]\t";
        savestream << "transVel[0]\t";
        savestream << "transVel[1]\t";
        savestream << "transVel[2]\t";
        savestream << "angVel[0]\t";
        savestream << "angVel[1]\t";
        savestream << "angVel[2]" << std::endl;        
    }
    
    savestream << step_id << "\t";
    savestream << sim_time << "\t";
    savestream << position[0] << "\t";
    savestream << position[1] << "\t";
    savestream << position[2] << "\t";
    savestream << quaternion[0] << "\t";
    savestream << quaternion[1] << "\t";
    savestream << quaternion[2] << "\t";
    savestream << quaternion[3] << "\t";
    savestream << transVel[0] << "\t";
    savestream << transVel[1] << "\t";
    savestream << transVel[2] << "\t";
    savestream << angVel[0] << "\t";
    savestream << angVel[1] << "\t";
    savestream << angVel[2] << std:: endl;
    savestream.close();
}*/
