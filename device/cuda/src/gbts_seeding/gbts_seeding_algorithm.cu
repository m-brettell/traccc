/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/utils.hpp"
#include "./kernels/GbtsNodesMakingKernels.cuh"
//#include "./kernels/GbtsGraphMaking.cuh"
//#include "./kernels/GbtsGraphProcessing.cuh"

#include "traccc/cuda/gbts_seeding/gbts_seeding_algorithm.hpp"

namespace traccc::cuda {
	
struct gbts_ctx {
	//counters
	unsigned int nSp{};
	unsigned int nEdges{};
	unsigned int nConnectedEdges{};
	unsigned int nSeeds{};
	//nEdges, nConnections, nConnectedEdges, .., nSeeds
	unsigned int* d_counters;	

	//NodeMaking
	unsigned int* d_layerCounts{};
	unsigned char* spacepointLayer{};
	//begin_idx for the surfaceToLayerMap or layerBin if one to one
	uint2* d_volumeToLayerMap{};	
	//surface_id, layerBin
	uint2* d_surfaceToLayerMap{};

	//x,y,z,cluster width in eta
	float4* d_reducedSP{};
	//output of layer binning
	float4* d_sp_params{};	
		
	//GraphMaking

	//GraphProccessing

	//output
}

gbts_seeding_algorithm::gbts_seeding_algorithm(const gbts_seedfinder_config& finder_config, traccc::memory_resource& mr, vecmem::copy& copy,stream& str, std::unique_ptr<const Logger> logger = getDummyLogger().clone())
	: messaging(std::move(logger), m_config(cfg), m_mr(mr), m_copy(copy) m_stream(str)) {}

output_type gbts_seeding::operator()(const edm::spacepoint_collection::const_view& spacepoints, const edm::measuremnt_collection::const_view& measurments) const {

	gbts_ctx ctx;

	//0. bin spacepoints by layer(disk) or any other maping supplied to the config.m_surfaceToLayerMap
	ctx.nSp = spacepoints.size();
	
	nThreads = 1024;
	nBlocks = std::ceil(ctx.nSp/nThreads);
	
	ctx.d_layerCounts = mr.main.do_allocate(mr, m_config.nLayers*sizeof(unsigned int));	
	ctx.d_spacepointsLayer = mr.main.do_allocate(mr, ctx.nSp*sizeof(unsigned char));	
	ctx.d_reducedSP = mr.main.do_allocate(mr, ctx.nSp*sizeof(float4));	
	
	ctx.d_volumeToLayerMap = mr.main.do_allocate(mr, sizeof(uint2)*m_config.nVolumes);	
	ctx.d_surfaceToLayerMap = mr.main.do_allocate(mr, sizeof(uint2)*m_config.nSurfaces);	
	
	m_copy.do_copy(sizeof(uint2)*m_config.nVolumes, m_config.volumeToLayerMap.data(), ctx.d_surfaceToLayerMap, cudaMemcpyHostToDevice);
	m_copy.do_copy(sizeof(uint2)*m_config.nSurface, m_config.surfaceToLayerMap.data(), ctx.d_surfaceToLayerMap, cudaMemcpyHostToDevice);
	
	count_sp_by_layer<<<nBlocks, nThreads, 0, m_stream>>>(spacepoints, measurments, ctx.d_volumeToLayerMap ,ctx.d_surfaceToLayerMap, ctx.d_layerCounts, ctx.d_spacepointsLayer, ctx.d_reducedSP, ctx.nSp);
	
	//prefix sum layerCounts
	m_copy.do_copy(m_config.nLayers*sizeof(unsigned int), ctx.d_layerCounts, layerPrefix, cudaMemcpyDeviceToHost);	
	std::unique_ptr<unsigned int[]> layerPrefix(m_config.nLayers);
	for(int layer = 1; later<m_config.nLayers; layer++) {
		layerPrefix[layer] += layerPrefix[layer-1]
	}
	m_copy.do_copy(m_config.nLayers*sizeof(unsigned int), layerPrefix, ctx.d_layerCounts, cudaMemcpyHostToDevice);	

	ctx.d_sp_params = mr.main.do_allocate(mr, ctx.nSp*sizeof(float4));	
	
	make_gbts_sp<<<nBlocks, nThreads, 0, m_stream>>>(ctx.d_sp_params, ctx.d_reducedSP, ctx.d_layerCount, ctx.d_spacepointLayer, ctx.nSp);	

	//1. histogram spacepoints by layer->eta->phi and convert to nodes phi,r,z,tau_min,tau_max

	//2. Find edges between spacepoint pairs

	//3. Link edges into graph
		
	//4. Prune unlinked edges from graph

	//5. Find longest segments with CCA

	//6. extract seeds, longest segment first
}
