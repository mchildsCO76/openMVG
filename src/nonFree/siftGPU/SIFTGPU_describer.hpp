// Copyright (c) 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef OPENMVG_PATENTED_SIFT_SIFTGPU_DESCRIBER_H
#define OPENMVG_PATENTED_SIFT_SIFTGPU_DESCRIBER_H

#include <cereal/cereal.hpp>
#include <iostream>
#include <numeric>
#include <string>

#include "nonFree/siftGPU/SiftGPU/src/SiftGPU/SiftGPU.h"
#include <GL/glut.h>
#include <GL/glu.h>


namespace openMVG {
namespace features {

// Bibliography:
// [1] R. Arandjelović, A. Zisserman.
// Three things everyone should know to improve object retrieval. CVPR2012.

inline void siftGPUDescToUChar(
  int n,
  std::vector<float> & descr,
  Descriptor<unsigned char,128> & descriptor,
  bool brootSift = false)
{
  if (brootSift)  {
    // rootsift = sqrt( sift / sum(sift) );
    const float sum = accumulate(descr.begin()+(n*128), descr.begin()+((n*128)+128), 0.0f);
    for (int k=0;k<128;++k){
      descriptor[k] = static_cast<unsigned char>(512.f * sqrt(descr[(n*128)+k]/sum));
    }
  }
  else{
    for (int k=0;k<128;++k){
		descriptor[k] = static_cast<unsigned char>(512.f*descr[(n*128)+k]);
	}
  }
}

struct SiftGPUParams
{
  SiftGPUParams(
    int first_octave = -2,
    int num_octaves = -2,
    int num_scales = -2,
    float edge_threshold = -2,
    float peak_threshold = -2,
    //int feat_level=2,
    //int feat_limit=7680,
    int feat_level=-2,
    int feat_limit=-2,
    bool adapt_dark = false,
    bool root_sift = true,
    bool use_cuda = true,
    int verbose_level = 0,
    std::string param_list = ""
  ):
    _first_octave(first_octave),
    _num_octaves(num_octaves),
    _num_scales(num_scales),
    _edge_threshold(edge_threshold),
    _peak_threshold(peak_threshold),
    _feat_level(feat_level),
    _feat_limit(feat_limit),
    _adapt_dark(adapt_dark),
    _root_sift(root_sift),
    _use_cuda(use_cuda),
    _verbose_level(verbose_level),
    _param_list(param_list) {}

  template<class Archive>
  void serialize( Archive & ar )
  {
    ar(
      cereal::make_nvp("first_octave", _first_octave),
      cereal::make_nvp("num_octaves",_num_octaves),
      cereal::make_nvp("num_scales",_num_scales),
      cereal::make_nvp("edge_threshold",_edge_threshold),
      cereal::make_nvp("peak_threshold",_peak_threshold),
      cereal::make_nvp("feat_level",_feat_level),
      cereal::make_nvp("feat_limit",_feat_limit),
      cereal::make_nvp("adapt_dark",_adapt_dark),
      cereal::make_nvp("root_sift",_root_sift),
      cereal::make_nvp("use_cuda",_use_cuda),
      cereal::make_nvp("verbose_level",_verbose_level),
      cereal::make_nvp("param_list",_param_list));
  }

  // Parameters
  // -2 means leave as default
  int _first_octave;      // Use original image, or perform an upscale if == -1
  int _num_octaves;       // Max octaves count
  int _num_scales;        // Scales per octave
  float _edge_threshold;  // Max ratio of Hessian eigenvalues
  float _peak_threshold;  // Min contrast
  int _feat_level;		  // Which levels to keep
  int _feat_limit;		  // Set a soft limit to number of detected features, provide -1 to disable
  bool _adapt_dark;
  bool _root_sift;        // see [1]
  bool _use_cuda;		  // Use CUDA (true) or GLSL (false)
  int _verbose_level;	  // Level of verbose
  std::string _param_list;	  // Set the parameters directly in the list
};

class SIFTGPU_Image_describer : public Image_describer
{
public:
  SIFTGPU_Image_describer(const SiftGPUParams & params = SiftGPUParams())
    :Image_describer(), _params(params){}

  ~SIFTGPU_Image_describer() {
	  delete sift_gpu;
  }

  bool Set_configuration_preset(EDESCRIBER_PRESET preset)
  {
	switch(preset)
    {
    case NORMAL_PRESET:
      _params._num_scales = 3 ;
    break;
    case HIGH_PRESET:
      _params._num_scales = 4;
    break;
    case ULTRA_PRESET:
      _params._num_scales = 7;
      _params._first_octave = -1;
    break;
    default:
      return false;
    }
    return true;
  }
  /**
  @brief SiftGPU detector has to be initialized before usage with the latest parameters set in _params
  */
  bool init(){
	  std::vector<std::string> params_vec;
    
    if(!_params._param_list.empty()){
		std::string param_string = _params._param_list;
		std::string delimiter=" ";
		std::string token;
		size_t pos = 0;
		while ((pos = param_string.find(delimiter)) != std::string::npos) {
			token = param_string.substr(0, pos);
			params_vec.push_back(token);
			param_string.erase(0, pos + delimiter.length());
		}
		params_vec.push_back(param_string);
	}
	else{
		if (_params._first_octave!=-2){
			params_vec.push_back(std::string("-fo"));
			params_vec.push_back(std::to_string(_params._first_octave));
		}
		if (_params._num_octaves!=-2){
			params_vec.push_back(std::string("-no"));
			params_vec.push_back(std::to_string(_params._num_octaves));
		}
		if (_params._num_scales!=-2){
			params_vec.push_back(std::string("-d"));
			params_vec.push_back(std::to_string(_params._num_scales));
		}
		if (_params._peak_threshold!=-2){
			params_vec.push_back(std::string("-t"));
			params_vec.push_back(std::to_string(_params._peak_threshold));
		}
		if (_params._edge_threshold!=-2){
			params_vec.push_back(std::string("-e"));
			params_vec.push_back(std::to_string(_params._edge_threshold));
		}
		if (_params._feat_level!=-2 && _params._feat_limit!=-2){
			params_vec.push_back(std::string("-tc").append(std::to_string(_params._feat_level)));
			params_vec.push_back(std::to_string(_params._feat_limit));
		}
		// Check if CUDA is supported
		#ifdef HAVE_SIFTGPU_CUDA
		if (_params._use_cuda){
			params_vec.push_back(std::string("-cuda"));
			params_vec.push_back("0");
			params_vec.push_back(std::string("-tc2"));
			params_vec.push_back("7680");
		}
		else if(_params._adapt_dark){
			params_vec.push_back(std::string("-da"));			
		}
		#else
		if(_params._adapt_dark){
			params_vec.push_back(std::string("-da"));			
		}
		#endif
		if (_params._verbose_level!=-2){
			params_vec.push_back(std::string("-v"));
			params_vec.push_back(std::to_string(_params._verbose_level));
		}

	}
    // Construct appropriate char array
    std::cout<<"--SiftGPU init params: ";
    char ** argv = new char*[params_vec.size()]();
    for (int p_i=0;p_i<params_vec.size();p_i++){
		argv[p_i] = &(params_vec[p_i][0]);
		std::cout<<params_vec[p_i]<<" ";
	}
    std::cout<<std::endl;
    
    if(params_vec.size()>0)
		sift_gpu->ParseParam(params_vec.size(), argv);
    
    //char * argv[] = {"-fo","-1","-t", "0.0002",  "-v", "1","-cuda ","0", "-m","1"};
    //sift_gpu->ParseParam(10, argv);
    
    if(sift_gpu->CreateContextGL() != SiftGPU::SIFTGPU_FULL_SUPPORTED) 
		return false;
		
    return true;
  }


  /**
  @brief Detect regions on the image and compute their attributes (description)
  @param image Image.
  @param regions The detected regions and attributes (the caller must delete the allocated data)
  @param mask 8-bit gray image for keypoint filtering (optional).
     Non-zero values depict the region of interest.
  */
  bool Describe(const image::Image<unsigned char>& image,
    std::unique_ptr<Regions> &regions,
    const image::Image<unsigned char> * mask = NULL)
  {
    const int w = image.Width(), h = image.Height();
    //Convert to float
    const image::Image<float> If(image.GetMat().cast<float>());
	//unsigned char* data = 
	sift_gpu->RunSIFT (w, h, If.data(), GL_LUMINANCE, GL_FLOAT);
	
	//get feature count
	int num = sift_gpu->GetFeatureNum();
	
	//allocate memory for readback
	vector<float> descriptors(128*num);
	vector<SiftGPU::SiftKeypoint> keys(num);
	
    Descriptor<unsigned char, 128> descriptor;
	
	//read back keypoints and normalized descritpros
	//specify NULL if you don’t need keypionts or descriptors
	sift_gpu->GetFeatureVector(&keys[0], &descriptors[0]);

    Allocate(regions);

    // Build alias to cached data
    SIFT_Regions * regionsCasted = dynamic_cast<SIFT_Regions*>(regions.get());
    regionsCasted->Features().reserve(num);
    regionsCasted->Descriptors().reserve(num);    

	for (int i = 0; i < num; i++) {

		// Feature masking
		if (mask)
		{
			const image::Image<unsigned char> & maskIma = *mask;
			if (maskIma(keys[i].y, keys[i].x) > 0)
				continue;
		}
		
    //Descriptor<unsigned char, 128> descriptor;
		// Create SIFT keypoint
		const SIOPointFeature fp(keys[i].x, keys[i].y,keys[i].s, keys[i].o);
		// Get 128 SIFT descriptor (it is normalized by 512)
		siftGPUDescToUChar(i, descriptors, descriptor, _params._root_sift);
		
		//Save keypoint and descriptor
		regionsCasted->Descriptors().push_back(descriptor);
		regionsCasted->Features().push_back(fp);
	}

    return true;
  };

  /// Allocate Regions type depending of the Image_describer
  void Allocate(std::unique_ptr<Regions> &regions) const
  {
    regions.reset( new SIFT_Regions );
  }

  template<class Archive>
  void serialize( Archive & ar )
  {
    ar(
     cereal::make_nvp("params", _params));
  }

private:
  SiftGPUParams _params;
  SiftGPU  *sift_gpu = new SiftGPU;
};

} // namespace features
} // namespace openMVG

#include <cereal/types/polymorphic.hpp>
#include <cereal/archives/json.hpp>
CEREAL_REGISTER_TYPE_WITH_NAME(openMVG::features::SIFTGPU_Image_describer, "SIFTGPU_Image_describer");

#endif // OPENMVG_PATENTED_SIFT_SIFT_DESCRIBER_H
