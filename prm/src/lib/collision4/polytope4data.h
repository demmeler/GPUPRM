#ifndef POLYTOPE4DATA_H
#define POLYTOPE4DATA_H

#include "cuda_head.h"
#include "polytope4.h"
#include "polytope.h"
#include "util.hpp"

#include "polytope4data.h"


namespace collision4{

#ifdef GPU_VERSION
    const int float4_storage_size=1000;
    __constant__ float4 float4_storage[float4_storage_size];
#endif

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


        //! data for pairs

        //! number of pairs (crs graph)
        //! only
        struct{
            int M;
            int *dsp;   //length N
            int *cnt;   //length N
            int *dest;  //length M
        }pairs;

        int *sys;

     public:
        //!get methods for kernel
        qualifier void get_polytope(polytope4& poly, int dof, int i) const;
        qualifier void get_polytope(polytope4& poly, int k) const;
        qualifier void get_collision_list(int k, const int* restrict &dest, int &num) const;
        qualifier int get_numsys(int dof) const{return numsys[dof];}

        //! init methods
        inline int build(const polytope4* polys, const int* sys, int N, int ndof_, const int* from_=0x0, const int* to_=0x0, int M_=0); //! arrays must be sortet w.r.t. sys
        inline int build(const polytope* polys, const int* sys, int N, int ndof_, const int* from_=0x0, const int* to_=0x0, int M_=0); //! arrays must be sortet w.r.t. sys

#ifdef GPU_VERSION
        friend int copy_host_to_device(polytope4data& devdata, const polytope4data& hostdata, bool withpairs_=false);
        friend int copy_host_to_device_ver2(polytope4data& devdata, const polytope4data& hostdata, bool withpairs_=false);
        friend int copy_device_to_host(polytope4data& hostdata, const polytope4data& devdata);
        //TODO: free memory
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

    qualifier void polytope4data::get_polytope(polytope4& poly, int k) const {
        poly.n=n[k];
        poly.m=m[k];
        int kn=dspn[k];
        poly.vertices=&vertices[kn];
        poly.dsp=&dsp[kn];
        poly.cnt=&cnt[kn];
        int km=dspm[k];
        poly.dest=&dest[km];
    }

    qualifier void polytope4data::get_collision_list(int k, const int* restrict &dest, int &num) const {
        num=pairs.cnt[k];
        dest=pairs.dest+pairs.dsp[k];
    }

    //! !! arrays must be sortet w.r.t. sys !!
    //! (from,to) must be sorted w.r.t. from!!
    inline int polytope4data::build(const polytope4* polys_, const int* sys_, int N_, int ndof_, const int* from_, const int* to_, int M_){
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
            if(k<dspsys[sys]){
                dspsys[sys]=k;
            }
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

        if(from_!=0x0 && to_!=0x0 && M_>0){
            pairs.dest=new int[M_];
            pairs.dsp=new int[N];
            pairs.cnt=new int[N]();
            pairs.M=M_;
            for(int i=0;i<M_;++i){
                pairs.dest[i]=to_[i];
                pairs.cnt[from_[i]]++;
            }
            int num=0;
            for(int k=0;k<N;++k){
                pairs.dsp[k]=num;
                num+=pairs.cnt[k];
            }
            check(num==M_);

            sys=new int[N];
            for(int i=0;i<N;++i){
                sys[i]=sys_[i];
            }
        }

        return 0; //Todo error handling?
    }


    //! !! arrays must be sortet w.r.t. sys !!
    //! (from,to) must be sorted w.r.t. from!!
    inline int polytope4data::build(const polytope* polys_, const int* sys_, int N_, int ndof_, const int* from_, const int* to_, int M_){
        ndof=ndof_;
        N=N_;
        n=new int[N];
        m=new int[N];
        dspn=new int[N];
        dspm=new int[N];
        dspsys=new int[ndof+1];
        for(int dof=0;dof<ndof+1;++dof){
            dspsys[dof]=N;
        }
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
            if(k<dspsys[sys]){
                dspsys[sys]=k;
            }
            numsys[sys]++;
        }

        vertices=new float4[sumn];
        dsp=new int[sumn];
        cnt=new int[sumn];
        dest=new int[summ];
        for(int k=0;k<N;++k){
            for(int i=0;i<n[k];++i){
                int l=dspn[k]+i;
                vertices[l].x=polys_[k].vertices[i].x;
                vertices[l].y=polys_[k].vertices[i].y;
                vertices[l].z=polys_[k].vertices[i].z;
                vertices[l].w=polys_[k].vertices[i].w;
                dsp[l]=polys_[k].dsp[i];
                cnt[l]=polys_[k].cnt[i];
            }
            for(int j=0;j<m[k];++j){
                int l=dspm[k]+j;
                dest[l]=polys_[k].dest[j];
            }
        }

        if(from_!=0x0 && to_!=0x0 && M_>0){
            pairs.dest=new int[M_];
            pairs.dsp=new int[N];
            pairs.cnt=new int[N]();
            pairs.M=M_;
            for(int i=0;i<M_;++i){
                pairs.dest[i]=to_[i];
                pairs.cnt[from_[i]]++;
            }
            int num=0;
            for(int k=0;k<N;++k){
                pairs.dsp[k]=num;
                num+=pairs.cnt[k];
            }
            check(num==M_);

            sys=new int[N];
            for(int i=0;i<N;++i){
                sys[i]=sys_[i];
            }
        }

        return 0; //Todo error handling?
    }

#ifdef GPU_VERSION
    inline int copy_host_to_device(polytope4data& devdata, const polytope4data& hostdata, bool withpairs){
        int N=devdata.N=hostdata.N;
        int sumn=devdata.sumn=hostdata.sumn;
        int summ=devdata.summ=hostdata.summ;
        int ndof=devdata.ndof=hostdata.ndof;

        //!alloc

        int res[28];
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


        if(withpairs){ //TODO: weglassen
            int M=devdata.pairs.M=hostdata.pairs.M;
            res[20]=cudaMalloc((void**)&(devdata.pairs.dsp), N * sizeof(int));
            res[21]=cudaMalloc((void**)&(devdata.pairs.cnt), N * sizeof(int));
            res[22]=cudaMalloc((void**)&(devdata.pairs.dest), M * sizeof(int));

            res[23]=cudaMemcpy((void*)devdata.pairs.dsp, (void*)hostdata.pairs.dsp, N * sizeof(int), cudaMemcpyHostToDevice);
            res[24]=cudaMemcpy((void*)devdata.pairs.cnt, (void*)hostdata.pairs.cnt, N * sizeof(int), cudaMemcpyHostToDevice);
            res[25]=cudaMemcpy((void*)devdata.pairs.dest, (void*)hostdata.pairs.dest, M * sizeof(int), cudaMemcpyHostToDevice);

            res[26]=cudaMalloc((void**)&(devdata.sys), N * sizeof(int));
            res[27]=cudaMemcpy((void*)devdata.sys, (void*)hostdata.sys, N * sizeof(int), cudaMemcpyHostToDevice);
        }


        for(int i=0;i<28;++i)if(res[i]!=0){printvar(i);printvar(res[i]); return 1000+i;}
        return 0;
    }

    inline int copy_device_to_host(polytope4data& hostdata, const polytope4data& devdata){
        msg("error: copy_host_to_device(polytope4data& hostdata, const polytope4data& devdata) not implemented");
        return -1;
    }


    inline int copy_host_to_device_ver2(polytope4data& devdata, const polytope4data& hostdata, bool withpairs){

        assert(hostdata.sumn<=float4_storage_size);

        int N=devdata.N=hostdata.N;
        int sumn=devdata.sumn=hostdata.sumn;
        int summ=devdata.summ=hostdata.summ;
        int ndof=devdata.ndof=hostdata.ndof;

        //!alloc

        cudaError_t res[28];
        //res[0]=cudaMalloc((void**)&(devdata.vertices.vertices), sumn * sizeof(float4));
        res[0]=cudaGetSymbolAddress((void**)&(devdata.vertices), float4_storage[0]);

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

        //res[10]=cudaMemcpy((void*)devdata.vertices, (void*)hostdata.vertices, sumn * sizeof(float4), cudaMemcpyHostToDevice);
        res[10]=cudaMemcpyToSymbol(float4_storage, (void*)hostdata.vertices, sumn * sizeof(float4), 0, cudaMemcpyHostToDevice);

        res[11]=cudaMemcpy((void*)devdata.dsp, (void*)hostdata.dsp, sumn * sizeof(int), cudaMemcpyHostToDevice);
        res[12]=cudaMemcpy((void*)devdata.cnt, (void*)hostdata.cnt, sumn * sizeof(int), cudaMemcpyHostToDevice);
        res[13]=cudaMemcpy((void*)devdata.dest, (void*)hostdata.dest, summ * sizeof(int), cudaMemcpyHostToDevice);

        res[14]=cudaMemcpy((void*)devdata.n, (void*)hostdata.n, N * sizeof(int), cudaMemcpyHostToDevice);
        res[15]=cudaMemcpy((void*)devdata.m, (void*)hostdata.m, N * sizeof(int), cudaMemcpyHostToDevice);
        res[16]=cudaMemcpy((void*)devdata.dspn, (void*)hostdata.dspn, N * sizeof(int), cudaMemcpyHostToDevice);
        res[17]=cudaMemcpy((void*)devdata.dspm, (void*)hostdata.dspm, N * sizeof(int), cudaMemcpyHostToDevice);

        res[18]=cudaMemcpy((void*)devdata.dspsys, (void*)hostdata.dspsys, (ndof+1) * sizeof(int), cudaMemcpyHostToDevice);
        res[19]=cudaMemcpy((void*)devdata.numsys, (void*)hostdata.numsys, (ndof+1) * sizeof(int), cudaMemcpyHostToDevice);


        if(withpairs){ //TODO: weglassen
            int M=devdata.pairs.M=hostdata.pairs.M;
            res[20]=cudaMalloc((void**)&(devdata.pairs.dsp), N * sizeof(int));
            res[21]=cudaMalloc((void**)&(devdata.pairs.cnt), N * sizeof(int));
            res[22]=cudaMalloc((void**)&(devdata.pairs.dest), M * sizeof(int));

            res[23]=cudaMemcpy((void*)devdata.pairs.dsp, (void*)hostdata.pairs.dsp, N * sizeof(int), cudaMemcpyHostToDevice);
            res[24]=cudaMemcpy((void*)devdata.pairs.cnt, (void*)hostdata.pairs.cnt, N * sizeof(int), cudaMemcpyHostToDevice);
            res[25]=cudaMemcpy((void*)devdata.pairs.dest, (void*)hostdata.pairs.dest, M * sizeof(int), cudaMemcpyHostToDevice);

            res[26]=cudaMalloc((void**)&(devdata.sys), N * sizeof(int));
            res[27]=cudaMemcpy((void*)devdata.sys, (void*)hostdata.sys, N * sizeof(int), cudaMemcpyHostToDevice);
        }


        for(int i=0;i<28;++i)if(res[i]!=0){printvar(i);printvar(res[i]);printvar(cudaGetErrorString(res[i])); return 1000+i;}
        return 0;
    }

#endif

}


#endif // POLYTOPE4DATA_H
