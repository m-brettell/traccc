/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/common.hpp"
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"

namespace traccc::device::GBTS {
//where to put these?
struct gbts_layerInfo {
	//int2 vector type
	int etaBin0{};
	int numBins{};
	//float2 vector type
	float minEta{};
	float delatEta{};
};

enum gbts_consts {
	
	//node binning
	unsigned short max_phi_bins       = 120,
	//graph building
	unsigned short node_buffer_length = 250,
	unsigned char  max_num_neighbours = 10,
	//graph processing
	unsigned char max_cca_iterations = 20,
	unsigned short shared_state_buffer_size = 578,
	//matrix access for kalman filter state
	unsigned char M3_0_0 = 0, M3_0_1 = 1, M3_0_2 = 2, M3_1_1 = 3, M3_1_2 = 4, M3_2_2 = 5,
		
	unsigned char M2_0_0 = 0, M2_0_1 = 1, M2_1_1 = 1,
}; 
}

namespace traccc {

struct gbts_seedfinder_config {
	gbts_seedfinder_config(const std::vector<std::pair<unsigned int, unsigned int>> binTables, std::vector<gbts_layerInfo> layerInfo,
              const std::vector<std::pair<uint64_t, unsigned int>> detrayBarcodeBinning, float minPt)
	
	//layer linking and geometry	
	std::vector<std::pair<unsigned char, unsigned char>> binTables{};
	std::vector<gbts_layerInfo> layerInfo{};
	unsigned int nLayers = 0;	
	
	std::vector<std::array<unsigned int, 2>> volumeToLayerMap;
	unsigned int nVolumes = 0;	

	std::vector<std::array<unsigned int, 2>> surfaceToLayerMap;
	unsigned int nSurfaces = 0;	

	//tuned for 900 MeV pT cut and scaled by input minPt	
	//edge making cuts
	float min_deltaPhi = 0.015;
	float dphi_coeff = 2.2e-4;
	float min_deltaPhi_low_dr = 0.002;
	float dphi_coeff_low_dr = 4.33e-4;	
	float minDeltaRadius = 2.0;
	float min_z0 = -600;
	float max_z0 = 600;
	float maxOuterRadius = 550;
	float cut_zMinU = 0;
	float cut_zMaxU = 0; //how to get ROI rzdr
	float maxKappa = 0.337;
	float low_Kappa_d0 = 0.02;
	float high_Kappa_d0 = 0.1;
	
	//edge matching cuts
	float cut_dphi_max = 0.012;
	float cut_dcurv_max = 0.001;
	float cut_tau_ratio_max = 0.01; 
	
	//seed extraction
};

//binTables contains pairs of linked layer-eta bins
//the layerInfo should be calculated from the barcodeBinning
//BarcodeBinning pair is detray barcode and bin index (corrisponding to the layers in layerInfo) ordered by volume 
//minPt in MeV
gbts_seedfinder_config::gbts_seedfinder_config(const std::vector<std::pair<unsigned int, unsigned char>> input_binTables, std::vector<gbts_layerInfo> input_layerInfo,
	  const std::vector<std::pair<uint64_t, unsigned int>> detrayBarcodeBinning, float minPt = 900) {
	
	//format layer-eta binning infomation
	binTables = input_binTables;
	layerInfo = input_layerInfo;
	unsigned int current_volume = 0;
	unsigned int current_layer = 0;
	bool layerChange = false;
	std::pair<uint64_t, unsigned int> surfaceLayerPair;
	for(unsigned int index = 0; index<detraySurfaceToLayerMap.size(); index++) {
		
		surfaceLayerPair = detraySurfaceToLayerMap[index];

		detray::geometry::barcode barcode(surfaceLayerPair.first);
		surfaceToLayerMap.push_back({barcode.surface_id, surfaceLayerPair.second});
		//is volume encompassed by a layer
		if(current_layer != surfaceLayerPair.second) layerChange = true;
		//create volume veiw on the surface map
		if(barcode.volume() != current_volume) {
			unsigned int bin = layerChange ? 1e6 : surfaceLayerPair.second;
			volumeToLayerMap.push_back({index, bin});
			
			layerChange = false;
		} 
	}
	//scale cuts
	float ptScale = 900/minPt;
	min_deltaPhi*=ptScale;
	dphi_coeff*=ptScale;
	min_deltaPhi_low_dr*=ptScale;
	dphi_coeff_low_dr*=ptScale;
	maxKappa*=ptScale;
	//contianers sizes
	nLayers   = layerInfo.size();
	nVolumes  = volumeToLayerMap.size();
	nSurfaces = surfaceToLayerMap.size();
}
} //namespace traccc
