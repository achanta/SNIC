//=================================================================================
//  snic.h
//
//
//  AUTORIGHTS
//  Copyright (C) 2019 Ecole Polytechnique Federale de Lausanne (EPFL), Switzerland.
//
//  Code released for research purposes only. For commercial purposes, please
//  contact the author at firstname.lastname@epfl.ch
//=================================================================================
/*Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met
 
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of EPFL nor the names of its contributors may
 be used to endorse or promote products derived from this software
 without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR EPFL BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */


typedef struct
{
    unsigned int i; // the x and y values packed into one
    unsigned int k; // the label
    double d;       // the distance
} NODE;

typedef struct
{
    NODE *nodes;
    int len; // number of elements present
    int size; // total capacity in terms of memory allocated
} HEAP;

void push (HEAP *h, const unsigned int ind, const unsigned int klab, const double dist);
// void pop (HEAP *h, NODE* pnode);
void pop(HEAP* h, unsigned int* ind, unsigned int* klab, double* dist);

void rgbtolab(double* rin, double* gin, double* bin, double sz, double* lvec, double* avec, double* bvec);
void FindSeeds(const int width, const int height, const int numk, int* kx, int* ky, int* outnumk);
void runSNIC(
            double**                      chans,
            const int                   nchans,
             const int                  width,
             const int                  height,
             int*                       labels,
             int*                       outnumk,
             const int                  innumk,
             const double               compactness);
void SNIC_main(double* img, const int width, const int height,
                const int nchannels, const int numSuperpixels, const double compactness,
                const int doRGBtoLAB, int* klabels, int* numlabels);

