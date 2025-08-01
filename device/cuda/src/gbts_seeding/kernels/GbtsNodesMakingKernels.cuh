/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

//cuda includes
#include <cuda.h>
#include <cuda_runtime.h>
#include <math_constants.h>
#include <vector_functions.h>

//Project includes
#include "traccc/cuda/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/edm/measurement.hpp"

// Detray include(s).
#include <detray/geometry/barcode.hpp>

namespace traccc::cuda::kernels {

__global__ void count_sp_by_layer(const traccc::edm::spacepoint_collection::const_view& spacpoints_view, const traccc::edm::measurment_collection_types::view& measurments_view, uint2 volumeToLayerMap, uint2 surfaceToLayerMap,
                                  unsigned int layerCounts, unsigned short spacepointsLayer, const unsigned int nSp, const unsigned int nSurfaces) {
	
	//copy map into shared?
	
	vecmem::device_vector<traccc::edm::measurement> measurements(measurements_view);
	vecmem::device_vector<traccc::edm::spacepoint> spacepoints(spacepoints_view);

	for(int spIdx = threadIdx.x + gridDim.x*blockDim.x; spIdx<nSp; spIdx += gridDim.x*blockDim.x) {
		//get the layer of the spacepoint
		traccc::edm::spacepoint spacepoint = spacepoints[spIdx];
		traccc::edm::measurment measurment = measurments[spacepoint.measurment_index_1()];
		
		detray::geometry::barcode barcode = measurment.surface_link;	
	
		unsigned int layerIdx;
		//some volume_ids map one to one with layer others need searching
		uint2 begin_or_bin = volumeToLayerMap[barcode.volume()];
		if(begin_or_bin.y == 1e6) {
			unsigned int surface_id = barcode.id();
			for(unsigned int surface = begin_or_bin.x; surface < nSurfaces ;surface++) {
				uint2 surfaceBinPair = surfaceToLayerMap[surface];
				if(surfaceBinPair.x == surface_id) { 
					layerIdx = surfaceBinPair.y; 
					break;
				}
			}
		}
		else layerIdx = begin_or_bin.y;
 
		//count and store x,y,z,cw info
		atomicAdd(&layerCounts[layerIdx], 1);
		spacepointsLayer[spIdx] = layerIdx;
		const traccc::point3 pos = spacepoint.global();
		reducedSP[spIdx] = make_float4(pos[0], pos[1], pos[2], measurment.diameter);
	} 	
}
 
//layerCount is prefix sumed on CPU inbetween count_sp_by_layer and this kerenel
__global__ void bin_sp_by_layer(float4 d_sp_params ,float4 reducedSP, unsigned int layerCount, unsigned short spacepointLayer, const unsigned int nSp) {
	
	for(int spIdx = threadIdx.x + gridDim.x*blockDim.x; spIdx<nSp; spIdx += gridDim.x*blockDim.x) {
		d_sp_params[atomicSub(&layerCount[spacepointLayer[spIdx]], 1)] = reducedSP[spIdx];
	}
}

__global__ void node_phi_binning_kernel(const float4* d_sp_params, int* d_node_phi_index, int nNodesPerBlock, int nNodes) {
    
static __constant__ int nBins = ;

    int begin_node = blockIdx.x * nNodesPerBlock;

    float inv_phiSliceWidth = 1/(2.0f * CUDART_PI_F/nBins);

    for(int idx = threadIdx.x + begin_node; idx < begin_node + nNodesPerBlock; idx += blockDim.x) {

       if (idx >= nNodes) continue;

       float4 sp = d_sp_params[idx];

       float Phi = atan2f(sp.y, sp.x);
   
       int phiIdx = (Phi + CUDART_PI_F)*inv_phiSliceWidth;
       
       if (phiIdx >= nBins) phiIdx %= nBins;
       else if (phiIdx < 0) {
            phiIdx += nBins;
            phiIdx %= nBins;
       }
       d_node_phi_index[idx] = phiIdx;
    }
}

__global__ void node_eta_binning_kernel(const float4* d_sp_params, const int4* d_layer_info, const float2* d_layer_geo,
                                        int* d_node_eta_index, int m_nLayers) {
    
    __shared__ int layer_begin;
    __shared__ int layer_end;
    __shared__ int num_eta_bins;
    __shared__ int bin0;
    __shared__ float min_eta;
    __shared__ float eta_bin_width;

    int layerIdx = blockIdx.x;

    if(threadIdx.x == 0) {
        int4 layerInfo     = d_layer_info[layerIdx];
        layer_begin        = layerInfo.x;
        layer_end          = layerInfo.y;
        num_eta_bins       = layerInfo.z;
        bin0               = layerInfo.w;
        float2 layerGeo    = d_layer_geo[layerIdx];
        min_eta            = layerGeo.x;
        eta_bin_width      = layerGeo.y;
    }

    __syncthreads();

    if(num_eta_bins == 1) {//all nodes are in the same bin
        for(int idx = threadIdx.x + layer_begin; idx < layer_end; idx += blockDim.x) {
            d_node_eta_index[idx] = bin0;
        }
    }
    else {//binIndex needs to be calculated first
        for(int idx = threadIdx.x + layer_begin; idx < layer_end; idx += blockDim.x) {

            float4 sp = d_sp_params[idx];

            float r   = sqrtf(powf(sp.x, 2) + powf(sp.y, 2));

            float t1   = sp.z/r;

            float eta = -logf(sqrtf(1 + powf(t1, 2)) - t1);

            int binIdx = (int)((eta - min_eta)/eta_bin_width);
  
            if(binIdx < 0) binIdx = 0;
            if(binIdx >= num_eta_bins) binIdx = num_eta_bins-1;
            
            d_node_eta_index[idx] = bin0 + binIdx;
        }
    }
}

__global__ void eta_phi_histo_kernel(const int* d_node_phi_index, const int* d_node_eta_index, unsigned int* d_eta_phi_histo, int nNodesPerBlock, int nNodes) {

    int begin_node = blockIdx.x * nNodesPerBlock;

    for(int idx = threadIdx.x + begin_node; idx < begin_node + nNodesPerBlock; idx += blockDim.x) {

       if (idx >= nNodes) continue;

       int eta_index = d_node_eta_index[idx];

       int histo_bin = d_node_phi_index[idx] + traccc::GBTS::num_phi_bins*eta_index;

       atomicAdd(&d_eta_phi_histo[histo_bin], 1);
    }

}

__global__ void eta_phi_counting_kernel(const unsigned int* d_histo, unsigned int* d_eta_node_counter, unsigned int* d_phi_cusums, int nBinsPerBlock, int maxEtaBin) {

    int eta_bin_start = nBinsPerBlock*blockIdx.x;

    int eta_bin_idx = eta_bin_start + threadIdx.x;
    
    if(eta_bin_idx >= maxEtaBin) return;

    int offset = traccc::GBTS::num_phi_bins*eta_bin_idx;
    
    int sum = 0;

    for(int phiIdx=0;phiIdx<traccc::GBTS::num_phi_bins;phiIdx++) {

        d_phi_cusums[offset + phiIdx] = sum;

        sum += d_histo[offset + phiIdx];
    }

    d_eta_node_counter[eta_bin_idx] = sum;
}

__global__ void eta_phi_prefix_sum_kernel(const unsigned int* d_histo, const unsigned int* d_eta_node_counter, unsigned int* d_phi_cusums, int nBinsPerBlock, int maxEtaBin) {

    int eta_bin_start = nBinsPerBlock*blockIdx.x;

    int eta_bin_idx = eta_bin_start + threadIdx.x;
    
    if(eta_bin_idx >= maxEtaBin) return;
    
    if(eta_bin_idx == 0) return;

    int offset = traccc::GBTS::num_phi_bins*eta_bin_idx;
    
    int val0 = d_eta_node_counter[eta_bin_idx-1];

    for(int phiIdx=0;phiIdx<traccc::GBTS::num_phi_bins;phiIdx++) {
        d_phi_cusums[offset + phiIdx] += val0;
    }
}

__global__ void node_sorting_kernel(const float4* d_sp_params, const int* d_node_eta_index, const int* d_node_phi_index, unsigned int* d_phi_cusums, 
                                    float* d_node_params, int* d_node_index, int nNodesPerBlock, int nNodes) {

    int begin_node = blockIdx.x * nNodesPerBlock;

    for(int idx = threadIdx.x + begin_node; idx < begin_node + nNodesPerBlock; idx += blockDim.x) {

       if (idx >= nNodes) continue;

       float4 sp = d_sp_params[idx];

       float Phi = atan2f(sp.y, sp.x);
       float r   = sqrtf(powf(sp.x, 2) + powf(sp.y, 2));
       float z = sp.z;

       float min_tau = -100.0;
       float max_tau = 100.0;

       if (sp.w > 0) {
            min_tau = 6.7*(sp.w - 0.2);//linear fit
            max_tau = 1.6 + 0.15/(sp.w + 0.2) + 6.1*(sp.w - 0.2);//linear fit + correction for short clusters
       }

       int eta_index = d_node_eta_index[idx];
       int histo_bin = d_node_phi_index[idx] + traccc::GBTS::num_phi_bins*eta_index;

       int pos = atomicAdd(&d_phi_cusums[histo_bin], 1);

       int o = 5*pos;

       d_node_params[o]   = min_tau;
       d_node_params[o+1] = max_tau;
       d_node_params[o+2] = Phi;
       d_node_params[o+3] = r;
       d_node_params[o+4] = z;
       d_node_index[pos] = idx;//keep the original index of the input spacepoint

    }                          
}

__global__ void minmax_rad_kernel(const int2* d_eta_bin_views, const float* d_node_params, float2* d_bin_rads, int nBinsPerBlock, int maxEtaBin) {

    int eta_bin_start = nBinsPerBlock*blockIdx.x;

    int eta_bin_idx = eta_bin_start + threadIdx.x;
    
    if(eta_bin_idx >= maxEtaBin) return;
    
    int2 view = d_eta_bin_views[eta_bin_idx];
    int node_start = view.x;
    int node_end = view.y;
    if (node_start == node_end) return;
    float min_r = 1e8;
    float max_r =-1e8;
    
    for(int idx = node_start; idx < node_end; idx++) {
        float r = d_node_params[5*idx + 3];
        if(r > max_r) max_r = r;
        if(r < min_r) min_r = r;
    }

    d_bin_rads[eta_bin_idx] = make_float2(min_r, max_r);

}
}
