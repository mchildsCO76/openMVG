

// Copyright (c) 2012, 2013, 2015 Pierre MOULON.

// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "openMVG/sfm/sfm.hpp"
#include "openMVG/image/image.hpp"

using namespace openMVG;
using namespace openMVG::features;
using namespace openMVG::cameras;
using namespace openMVG::geometry;
using namespace openMVG::image;
using namespace openMVG::sfm;

#include "third_party/cmdLine/cmdLine.h"
#include "third_party/stlplus3/filesystemSimplified/file_system.hpp"

#include <fstream>

void getLandmarksPerView(
  const SfM_Data & sfm_data,
  Hash_Map< IndexT, std::vector<Vec3> > & landmarks_per_view
  )
{
  for (Landmarks::const_iterator iterTracks = sfm_data.GetLandmarks().begin();
    iterTracks != sfm_data.GetLandmarks().end();
    ++iterTracks
  )
  {
    const Observations & obs = iterTracks->second.obs;
    for (Observations::const_iterator itObs = obs.begin();
      itObs != obs.end(); ++itObs)
    {
      // Add Z value of landmark
      landmarks_per_view[itObs->first].push_back(iterTracks->second.X);
    }
  }

}

bool exportCAMSToMICMAC(
  const SfM_Data & sfm_data,
  const std::string & sOutDir,
  const bool bExportImages)
{
  std::ostringstream os;
  
  // Create the Orientation folder
  const std::string orientation_str = "Ori-OpenMVG";
  const std::string & sOutCalibFolder = stlplus::folder_append_separator(sOutDir) + orientation_str;
  if (!stlplus::is_folder(sOutCalibFolder)){
    stlplus::folder_create(sOutCalibFolder);
  }

  // Check of folders are correctly created
  bool bOk = (stlplus::is_folder(sOutCalibFolder));

  if(bOk){
    // Create MicMac-LocalChantierDescripteur.xml file
    const std::string & sLocalChantierDescripteur = stlplus::folder_append_separator(sOutDir)+"MicMac-LocalChantierDescripteur.xml";
    std::ofstream file_MM_lcd_xml(sLocalChantierDescripteur.c_str());
    file_MM_lcd_xml << "<Global>" << os.widen('\n')
     << os.widen('\t')<<"<ChantierDescripteur >"<< os.widen('\n');


    // Compute Landmarks for each view
    Hash_Map< IndexT, std::vector<Vec3> > landmarks_per_view;
    getLandmarksPerView(sfm_data, landmarks_per_view);


    std::vector<Views::const_iterator> valid_views;
    // Export intrinsic parameters for each camera
    for(Views::const_iterator iter = sfm_data.GetViews().begin();
    iter != sfm_data.GetViews().end(); ++iter)
    {
      const View * view = iter->second.get();
      const IndexT id_view = view->id_view;

      if (!sfm_data.IsPoseAndIntrinsicDefined(view))
        continue;

      // Valid view, we can ask a pose & intrinsic data
      const Pose3 pose = sfm_data.GetPoseOrDie(view);
      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);
      const IntrinsicBase * cam = iterIntrinsic->second.get();

      if (!cameras::isPinhole(cam->getType()))
        continue;

      const Pinhole_Intrinsic * pinhole_cam = static_cast<const Pinhole_Intrinsic *>(cam);

      // We compute the 'AltiSol' an the 'Depth' param. Suboptimal for now.
      double meanDepth = 0., meanZ = 0.;

      // Find all landmarks it sees
      if (landmarks_per_view.find(id_view)!=landmarks_per_view.end())
      {
        // Iterate and get average Z value of the landmarks
        std::vector<Vec3> & landmarks_view = landmarks_per_view.at(id_view);
        int n_landmarks=0;
        for(std::vector<Vec3>::iterator it = landmarks_view.begin(); it!= landmarks_view.end(); ++it){
          n_landmarks++;
          // Use running average to avoid possibility of overflow
          meanZ = (meanZ * (n_landmarks-1)/(n_landmarks*1.0)) + ((*it)(2)/(n_landmarks*1.0));
          meanDepth = (meanDepth * (n_landmarks-1)/(n_landmarks*1.0)) + (pose.depth((*it))/(n_landmarks*1.0));
        }
      }
      valid_views.push_back(iter);

      // Extrinsic
      const Vec3 c = pose.center(); // both have center in world c.s w_t
      const Mat3 R = pose.rotation().transpose(); // MICMAC has w_R_c while openMVG c_R_w

      // Intrinsics ID
      IndexT intrinsic_id = view->id_intrinsic;
      // Intrinsic
      const double f = pinhole_cam->focal();
      const int focal_int = (int)(f+0.05)*10;

      const Vec2 pp = pinhole_cam->principal_point();

      // Image size in px
      const int w = pinhole_cam->w();
      const int h = pinhole_cam->h();

      // Create orientation file
      std::ofstream outfile( stlplus::create_filespec(
           sOutCalibFolder, std::string("Orientation-") + stlplus::filename_part(view->s_Img_path), "xml" ).c_str() );
      outfile << setprecision(8) << fixed ;

      // See ParamChantierPhotogram.xml in MicMac distrib for full specs. The doc is also useful !
      outfile
        << "<?xml version=\"1.0\" ?>\n"
        << "<ExportAPERO>\n"
        << "    <OrientationConique>\n"
        << "        <OrIntImaM2C>\n"
        << "            <I00>0 0</I00>\n"
        << "            <V10>1 0</V10>\n"
        << "            <V01>0 1</V01>\n"
        << "        </OrIntImaM2C>\n"
        << "        <TypeProj>eProjStenope</TypeProj>\n"
        << "        <ZoneUtileInPixel>true</ZoneUtileInPixel>\n"
        << "        <Interne>\n"
        << "            <KnownConv>eConvApero_DistM2C</KnownConv>\n"
        << "            <PP>" << pp(0) << " " << pp(1) <<"</PP>\n"
        << "            <F>" << f << "</F>\n"
        << "            <SzIm>" << w << " " << h << "</SzIm>\n"
        << "            <CalibDistortion>\n"
        << "                <ModRad>\n"
        << "                    <CDist>" << pp(0) << " " <<  pp(1) << "</CDist>\n"
        << "                </ModRad>\n"
        << "            </CalibDistortion>\n"
        << "        </Interne>\n"
        << "        <Externe>\n"
        << "            <AltiSol>" << meanZ << "</AltiSol>\n"
        << "            <Profondeur>" << meanDepth << "</Profondeur>\n"
        << "            <KnownConv>eConvApero_DistM2C</KnownConv>\n"
        << "            <Centre>" << c(0) << " " << c(1) << " " << c(2) << "</Centre>\n"
        << "            <IncCentre>1 1 1</IncCentre>\n"
        << "            <ParamRotation>\n"
        << "                <CodageMatr>\n"
        << "                    <L1>" << R(0,0) << " " << R(0,1) << " " << R(0,2) <<  "</L1>\n"
        << "                    <L2>" << R(1,0) << " " << R(1,1) << " " << R(1,2) << "</L2>\n"
        << "                    <L3>" << R(2,0) << " " << R(2,1) << " " << R(2,2) << "</L3>\n"
        << "                </CodageMatr>\n"
        << "            </ParamRotation>\n"
        << "        </Externe>\n"
        << "        <ConvOri>\n"
        << "            <KnownConv>eConvApero_DistM2C</KnownConv>\n"
        << "        </ConvOri>\n"
        << "    </OrientationConique>\n"
        << "</ExportAPERO>\n";
      outfile.close();

        // Copy images to output folder
      if(bExportImages){
        Image<RGBColor> image, image_ud;
        const std::string srcImage = stlplus::create_filespec(sfm_data.s_root_path,view->s_Img_path);
        const std::string dstImage = stlplus::create_filespec(sOutDir,stlplus::filename_part(view->s_Img_path) ,"");
        if (cam->have_disto())
        {
        // undistort the image and save it
          if (ReadImage( srcImage.c_str(), &image))
          {
            UndistortImage(image, cam, image_ud, BLACK);
            bOk &= WriteImage(dstImage.c_str(), image_ud);
          }
        }
        else // (no distortion)
        {
          // copy the image since there is no distortion
          stlplus::file_copy(srcImage, dstImage);
        }
      }

    }


    file_MM_lcd_xml << os.widen('\t') << "<KeyedNamesAssociations>" << os.widen('\n');
    for(int i=0;i<valid_views.size();i++){
      const View * view = valid_views.at(i)->second.get();
      std::string filename = stlplus::filename_part(view->s_Img_path);
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << "<Calcs>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << "<Arrite> 1 1 </Arrite>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << "<Direct>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << os.widen('\t') << "<PatternTransform>"<<filename<< "</PatternTransform>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << os.widen('\t') << "<CalcName>"<<"OpenMVG_cam_"<<view->id_intrinsic<< "</CalcName>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << "</Direct>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << "</Calcs>" << os.widen('\n');
    }
    file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << "<Key>   NKS-Assoc-STD-CAM  </Key>" << os.widen('\n');
    file_MM_lcd_xml << os.widen('\t') << "</KeyedNamesAssociations>" << os.widen('\n');

    file_MM_lcd_xml << os.widen('\t') << "<KeyedNamesAssociations>" << os.widen('\n');
    for(int i=0;i<valid_views.size();i++){
      const View * view = valid_views.at(i)->second.get();
      std::string filename = stlplus::filename_part(view->s_Img_path);

      Intrinsics::const_iterator iterIntrinsic = sfm_data.GetIntrinsics().find(view->id_intrinsic);
      const double f = static_cast<const Pinhole_Intrinsic *>((sfm_data.GetIntrinsics().find(view->id_intrinsic))->second.get())->focal();
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << "<Calcs>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << "<Arrite> 1 1 </Arrite>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << "<Direct>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << os.widen('\t') << "<PatternTransform>"<<filename<< "</PatternTransform>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << os.widen('\t') << "<CalcName>"<<f<< "</CalcName>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << os.widen('\t') << "</Direct>" << os.widen('\n');
      file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << "</Calcs>" << os.widen('\n');
    }
    file_MM_lcd_xml << os.widen('\t') << os.widen('\t') << "<Key>   NKS-Assoc-STD-FOC  </Key>" << os.widen('\n');
    file_MM_lcd_xml << os.widen('\t') << "</KeyedNamesAssociations>" << os.widen('\n');

    file_MM_lcd_xml << os.widen('\t') << "</ChantierDescripteur>" << os.widen('\n');
    file_MM_lcd_xml << os.widen('\t') << "</Global>" << os.widen('\n');
    file_MM_lcd_xml.close();

  }
  else{
    return false;
  }
  return true;

}

bool exportPOINTSToMICMAC(
  const SfM_Data & sfm_data,
  const std::string & sOutDir,
  const bool bExportImages)
{
  std::ostringstream os;

  // Create the Orientation folder
  const std::string homol_str = "Homol";
  const std::string & sOutHomolFolder = stlplus::folder_append_separator(sOutDir) + homol_str;
  if (!stlplus::is_folder(sOutHomolFolder)){
    stlplus::folder_create(sOutHomolFolder);
  }

  // Check of folders are correctly created
  bool bOk = (stlplus::is_folder(sOutHomolFolder));

  if(bOk){
    Hash_Map<IndexT, Hash_Map<IndexT, std::vector<std::pair<Observation,Observation> > > > obs_per_view_pair;

      std::cout << "\n" << "Landmark generation" << std::endl;
    C_Progress_display pairs_progress_bar( sfm_data.GetLandmarks().size() );
    for (Landmarks::const_iterator iterTracks = sfm_data.GetLandmarks().begin();
    iterTracks != sfm_data.GetLandmarks().end();
    ++iterTracks
    )
    {
      const Observations & obs = iterTracks->second.obs;
      //Get through observations and link all of them to all
      for (Observations::const_iterator itObs_A = obs.begin(); itObs_A != obs.end(); ++itObs_A)
      {
        IndexT view_A = itObs_A->first;
        Observation obs_A = itObs_A->second;

        Observations::const_iterator itObs_B =  std::next(itObs_A, 1);
        for (; itObs_B != obs.end(); ++itObs_B)
        {
          IndexT view_B = itObs_B->first;
          Observation obs_B = itObs_B->second;

          if(view_B < view_A){
            IndexT tmp_i = view_A;
            Observation tmp_o = obs_A;
            view_A = view_B;
            obs_A = obs_B;
            view_B = tmp_i;
            obs_B = tmp_o;
          }

          if(obs_per_view_pair.find(view_A)!=obs_per_view_pair.end()){
            Hash_Map<IndexT, std::vector<std::pair<Observation,Observation> > > &obs_per_view = obs_per_view_pair.at(view_A);
            if(obs_per_view.find(view_B)!=obs_per_view.end()){
              std::vector<std::pair<Observation,Observation> > & obs_vector =  obs_per_view.at(view_B);
              obs_vector.push_back(std::pair<Observation,Observation>(obs_A,obs_B));
            }
            else{
              std::vector<std::pair<Observation,Observation> > obs_vector;
              obs_vector.push_back(std::pair<Observation,Observation>(obs_A,obs_B));
              obs_per_view.insert(std::pair<IndexT, std::vector<std::pair<Observation,Observation> > >(view_B,obs_vector));
            }
          }
          else{
            Hash_Map<IndexT, std::vector<std::pair<Observation,Observation> > > obs_per_view;
            std::vector<std::pair<Observation,Observation> > obs_vector;
            obs_vector.push_back(std::pair<Observation,Observation>(obs_A,obs_B));
            obs_per_view.insert(std::pair<IndexT, std::vector<std::pair<Observation,Observation> > >(view_B,obs_vector));
            obs_per_view_pair.insert(std::pair<IndexT, Hash_Map<IndexT, std::vector<std::pair<Observation,Observation> > > >(view_A,obs_per_view));
          }

        }
      }
      ++pairs_progress_bar;
    }

    std::cout << "\n Export landmarks" << std::endl;
    pairs_progress_bar.restart( obs_per_view_pair.size() );

    for(Hash_Map<IndexT, Hash_Map<IndexT, std::vector<std::pair<Observation,Observation> > > >::iterator view_A_iter = obs_per_view_pair.begin();view_A_iter!=obs_per_view_pair.end();++view_A_iter){

      IndexT view_A_id = view_A_iter->first;
      View * view_A = sfm_data.GetViews().at(view_A_id).get();

      Hash_Map<IndexT, std::vector<std::pair<Observation,Observation> > > obs_per_view_A = view_A_iter->second;

      for(Hash_Map<IndexT, std::vector<std::pair<Observation,Observation> > >::iterator view_B_iter = obs_per_view_A.begin();view_B_iter!=obs_per_view_A.end();++view_B_iter){

        IndexT view_B_id = view_B_iter->first;
        View * view_B = sfm_data.GetViews().at(view_B_id).get();

        std::vector<std::pair<Observation,Observation> >  obs_vector =  view_B_iter->second;

        const std::string filename_A = stlplus::filename_part(view_A->s_Img_path);
        const std::string filename_B = stlplus::filename_part(view_B->s_Img_path);

        const std::string & sImgA_folder_path = stlplus::folder_append_separator(sOutHomolFolder)+"Pastis"+ filename_A;
        const std::string & sImgB_folder_path = stlplus::folder_append_separator(sOutHomolFolder)+"Pastis"+ filename_B;
        const std::string & sImgB_A_file_path = stlplus::folder_append_separator(sImgB_folder_path)+""+ filename_A+".dat";
        const std::string & sImgA_B_file_path = stlplus::folder_append_separator(sImgA_folder_path)+""+ filename_B+".dat";

        // Create export folders
        if (!stlplus::is_folder(sImgA_folder_path)){
          stlplus::folder_create(sImgA_folder_path);
        }
        if (!stlplus::is_folder(sImgA_folder_path)){
          stlplus::folder_create(sImgA_folder_path);
        }

        int n_obs = obs_vector.size();
        int num_2 = 2;
        double num_1 = 1;


        // Write header to each file
        std::ofstream file_img_A_B(sImgA_B_file_path.c_str(),std::ofstream::binary);
        file_img_A_B.write(reinterpret_cast<const char *>(&num_2),sizeof(num_2));
        file_img_A_B.write(reinterpret_cast<const char *>(&n_obs),sizeof(n_obs));
        std::ofstream file_img_B_A(sImgB_A_file_path.c_str(),std::ofstream::binary);
        file_img_B_A.write(reinterpret_cast<const char *>(&num_2),sizeof(num_2));
        file_img_B_A.write(reinterpret_cast<const char *>(&n_obs),sizeof(n_obs));

        for (const std::pair< Observation, Observation >& point_pair : obs_vector)
        {
          // Write the feature pair
          file_img_A_B.write(reinterpret_cast<const char *>(&num_2),sizeof(num_2));
          file_img_A_B.write(reinterpret_cast<const char *>(&num_1),sizeof(num_1));
          file_img_B_A.write(reinterpret_cast<const char *>(&num_2),sizeof(num_2));
          file_img_B_A.write(reinterpret_cast<const char *>(&num_1),sizeof(num_1));

          const double im_A_x =  point_pair.first.x(0);
          const double im_A_y =  point_pair.first.x(1);
          const double im_B_x =  point_pair.second.x(0);
          const double im_B_y =  point_pair.second.x(1);

          file_img_A_B.write(reinterpret_cast<const char *>(&im_A_x),sizeof(im_A_x));
          file_img_A_B.write(reinterpret_cast<const char *>(&im_A_y),sizeof(im_A_y));
          file_img_A_B.write(reinterpret_cast<const char *>(&im_B_x),sizeof(im_B_x));
          file_img_A_B.write(reinterpret_cast<const char *>(&im_B_y),sizeof(im_B_y));

          file_img_B_A.write(reinterpret_cast<const char *>(&im_B_x),sizeof(im_B_x));
          file_img_B_A.write(reinterpret_cast<const char *>(&im_B_y),sizeof(im_B_y));
          file_img_B_A.write(reinterpret_cast<const char *>(&im_A_x),sizeof(im_A_x));
          file_img_B_A.write(reinterpret_cast<const char *>(&im_A_y),sizeof(im_A_y));
        }
        file_img_A_B.close();
        file_img_B_A.close();

      }
      ++pairs_progress_bar;
    }
  }

}

int main(int argc, char **argv)
{
  CmdLine cmd;

  std::string sSfM_Data_Filename;
  std::string sOutDir = "";
  bool bExportImages = true;
  bool bExportFeatures = true;

  cmd.add( make_option('i', sSfM_Data_Filename, "sfmdata") );
  cmd.add( make_option('o', sOutDir, "outdir") );
  cmd.add( make_option('f', bExportFeatures, "exportFeatures") );
  cmd.add( make_option('e', bExportImages, "exportImages") );

  try {
    if (argc == 1) throw std::string("Invalid command line parameter.");
    cmd.process(argc, argv);
  } catch(const std::string& s) {
    std::cerr << "Usage: " << argv[0] << '\n'
    << "[-i|--sfmdata] filename, the SfM_Data (after SfM) file to convert\n"
    << "[-o|--outdir] path.\n"
    << "[-e|--exportImages] \n"
    << "\t 0-> do NOT export images\n"
    << "\t 1-> export images to output folder (default).\n"
    << "[-f|--exportFeatures] \n"
    << "\t 0-> do NOT export features\n"
    << "\t 1-> export features (default).\n"
    << std::endl;

    std::cerr << s << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << " You called : " <<std::endl
      << argv[0] << std::endl
      << "--sfmdata " << sSfM_Data_Filename << std::endl
      << "--outdir " << sOutDir << std::endl
      << "--exportFeatures " << bExportFeatures << std::endl
      << "--exportImages " << bExportImages << std::endl;

  // Create output dir
  if (!stlplus::folder_exists(sOutDir))
    stlplus::folder_create( sOutDir );

  // Read the SfM scene
  SfM_Data sfm_data;
  if (!Load(sfm_data, sSfM_Data_Filename, ESfM_Data(ALL))) {
    std::cerr << std::endl
      << "The input SfM_Data file \""<< sSfM_Data_Filename << "\" cannot be read." << std::endl;
    return EXIT_FAILURE;
  }


  // Export camera orientations
  exportCAMSToMICMAC(sfm_data,sOutDir,bExportImages);
  // Export features
  if(bExportFeatures)
    exportPOINTSToMICMAC(sfm_data,sOutDir,bExportImages);

  return EXIT_SUCCESS;
}