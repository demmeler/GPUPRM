#ifndef POLYTOPE4DATA_H
#define POLYTOPE4DATA_H

#include "cuda_head.h"
#include "polytope4.h"
#include "util.hpp"

namespace collision4{

    class polytope4data{
    public:
        //!raw data of all polytopes

        float4 *vertices; //length sumn
        int *dsp;         //    "
        int *cnt;         //    "
        int *n;           //length N

        int *dest;        //length summ
        int *m;           //length N

        //!organization

        //! sum over all n
        int sumn;
        //! sum over all m
        int summ;

        //!number of polytopes
        int N;
        //!displacement arrays for sum(n) arrays
        int *dspn;        //length N
        //!and sum(m) array
        int *dspm;        //length N

        //!assignment to dof system
        int *dspsys;      //length ndof+1
        int *numsys;      //length ndof+1

        //!number degrees of freedom
        int ndof;

     public:
        //!get methods for kernel
        qualifier void get_polytope(polytope4& poly, int dof, int k) const;
        qualifier int get_numsys(int dof) const{return numsys[dof];}

        //! init methods
        int build(polytope4* polys, int* sys, int N, int ndof_); //! arrays must be sortet w.r.t. dof
#ifdef CUDA_IMPLEMENTATION
        friend int copy_host_to_device(polytope4data& devdata, const polytope4data& hostdata);
        friend int copy_device_to_host(polytope4data& hostdata, const polytope4data& devdata);
#endif

    };



    qualifier void polytope4data::get_polytope(polytope4& poly, int dof, int i) const {
        int k=dspsys[dof]+i;
        poly.n=n[k];
        poly.m=m[k];
        int kn=dspn[k];
        poly.vertices=&vertices[kn];
        poly.dsp=&dsp[kn];
        poly.cnt=&cnt[kn];
        int km=dspm[k];
        poly.dest=&dest[km];
    }


    //! arrays must be sortet w.r.t. dof
    int polytope4data::build(polytope4* polys_, int* sys_, int N_, int ndof_){
        ndof=ndof_;
        N=N_;
        n=new int[N];
        m=new int[N];
        dspn=new int[N];
        dspm=new int[N];
        dspsys=new int[ndof+1];
        for(int dof=0;dof<ndof+1;++dof){dspsys[dof]=N;}
        numsys=new int[ndof+1]();
        sumn=summ=0;
        for(int k=0;k<N;++k){
            n[k]=polys_[k].n;
            m[k]=polys_[k].m;
            dspn[k]=sumn;
            dspm[k]=summ;
            sumn+=n[k];
            summ+=m[k];
            int sys=sys_[k];
            dspsys[sys]=k;
            numsys[sys]++;
        }

        vertices=new float4[sumn];
        dsp=new int[sumn];
        cnt=new int[sumn];
        dest=new int[summ];
        for(int k=0;k<N;++k){
            for(int i=0;i<n[k];++i){
                int l=dspn[k]+i;
                vertices[l]=polys_[k].vertices[i];
                dsp[l]=polys_[k].dsp[i];
                cnt[l]=polys_[k].cnt[i];
            }
            for(int j=0;j<m[k];++j){
                int l=dspm[k]+j;
                dest[l]=polys_[k].dest[j];
            }
        }
        return 0; //Todo error handling?
    }

#ifdef CUDA_IMPLEMENTATION
    int copy_host_to_device(polytope4data& devdata, const polytope4data& hostdata){
        int N=devdata.N=hostdata.N;
        int sumn=devdata.sumn=hostdata.sumn;
        int summ=devdata.summ=hostdata.summ;
        int ndof=devdata.ndof=hostdata.ndof;

        //!alloc

        int res[20];
        res[0]=cudaMalloc((void**)&(devdata.vertices), sumn * sizeof(float4));
        res[1]=cudaMalloc((void**)&(devdata.dsp), sumn * sizeof(int));
        res[2]=cudaMalloc((void**)&(devdata.cnt), sumn * sizeof(int));
        res[3]=cudaMalloc((void**)&(devdata.dest), summ * sizeof(int));

        res[4]=cudaMalloc((void**)&(devdata.n), N * sizeof(int));
        res[5]=cudaMalloc((void**)&(devdata.m), N * sizeof(int));
        res[6]=cudaMalloc((void**)&(devdata.dspn), N * sizeof(int));
        res[7]=cudaMalloc((void**)&(devdata.dspm), N * sizeof(int));

        res[8]=cudaMalloc((void**)&(devdata.dspsys),(ndof+1) * sizeof(int));
        res[9]=cudaMalloc((void**)&(devdata.numsys),(ndof+1) * sizeof(int));

        //! copy

        res[10]=cudaMemcpy((void*)devdata.vertices, (void*)hostdata.vertices, sumn * sizeof(float4), cudaMemcpyHostToDevice);
        res[11]=cudaMemcpy((void*)devdata.dsp, (void*)hostdata.dsp, sumn * sizeof(int), cudaMemcpyHostToDevice);
        res[12]=cudaMemcpy((void*)devdata.cnt, (void*)hostdata.cnt, sumn * sizeof(int), cudaMemcpyHostToDevice);
        res[13]=cudaMemcpy((void*)devdata.dest, (void*)hostdata.dest, summ * sizeof(int), cudaMemcpyHostToDevice);

        res[14]=cudaMemcpy((void*)devdata.n, (void*)hostdata.n, N * sizeof(int), cudaMemcpyHostToDevice);
        res[15]=cudaMemcpy((void*)devdata.m, (void*)hostdata.m, N * sizeof(int), cudaMemcpyHostToDevice);
        res[16]=cudaMemcpy((void*)devdata.dspn, (void*)hostdata.dspn, N * sizeof(int), cudaMemcpyHostToDevice);
        res[17]=cudaMemcpy((void*)devdata.dspm, (void*)hostdata.dspm, N * sizeof(int), cudaMemcpyHostToDevice);

        res[18]=cudaMemcpy((void*)devdata.dspsys, (void*)hostdata.dspsys, (ndof+1) * sizeof(int), cudaMemcpyHostToDevice);
        res[19]=cudaMemcpy((void*)devdata.numsys, (void*)hostdata.numsys, (ndof+1) * sizeof(int), cudaMemcpyHostToDevice);

        for(int i=0;i<20;++i)if(res[i]!=0)return res[i];
        return 0;
    }

    int copy_device_to_host(polytope4data& hostdata, const polytope4data& devdata){
        msg("error: copy_host_to_device(polytope4data& hostdata, const polytope4data& devdata) not implemented");
        return -1;
    }
#endif

}


#endif // POLYTOPE4DATA_H
