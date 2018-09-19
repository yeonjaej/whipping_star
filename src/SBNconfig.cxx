#include "SBNconfig.h"
using namespace sbn;

//standard constructor given an .xml
SBNconfig::SBNconfig(std::string whichxml, bool isverbose): xmlname(whichxml) {
    otag = "SBNconfig::SBNconfig\t|| ";

    is_verbose = isverbose;
    has_oscillation_patterns = false;

    //max subchannels 100?
    subchannel_names.resize(100);
    subchannel_bool.resize(100);
    subchannel_osc_patterns.resize(100);
    char *end;

    //Setup TiXml documents
    TiXmlDocument doc( whichxml.c_str() );
    bool loadOkay = doc.LoadFile();
    if(loadOkay){
        if(is_verbose)	std::cout<<otag<<"Loaded XML configuration file: "<<whichxml<<std::endl;
    }else{
        std::cout<<otag<<"ERROR: Failed to load XML configuration file:  "<<whichxml<<std::endl;
        exit(EXIT_FAILURE);
    }
    TiXmlHandle hDoc(&doc);


    // we have Modes, Detectors, Channels
    TiXmlElement *pMode, *pDet, *pChan, *pCov, *pMC, *pData;


    //Grab the first element. Note very little error checking here! make sure they exist.
    pMode = doc.FirstChildElement("mode");
    pDet =  doc.FirstChildElement("detector");
    pChan = doc.FirstChildElement("channel");
    pCov  = doc.FirstChildElement("covariance");
    pMC   = doc.FirstChildElement("MultisimFile");
    pData   = doc.FirstChildElement("data");

    if(!pMode){
        std::cout<<otag<<"ERROR: Need at least 1 mode defined in xml./n";
    }else{
        while(pMode){
            // What modes are we running in (e.g nu, nu bar, horn current=XXvolts....) Can have as many as we want
            mode_names.push_back(pMode->Attribute("name"));
            mode_bool.push_back(strtod(pMode->Attribute("use"),&end));

            pMode = pMode->NextSiblingElement("mode");
            if(is_verbose)	std::cout<<"SBNconfig::SBnconfig\t|| loading mode: "<<mode_names.back()<<" with use_bool "<<mode_bool.back()<<std::endl;

        }
    }

    // How many detectors do we want!
    if(!pDet){
        std::cout<<otag<<"ERROR: Need at least 1 detector defined in xml./n";
    }else{

        while(pDet){
            //std::cout<<"Detector: "<<pDet->Attribute("name")<<" "<<pDet->Attribute("use")<<std::endl;
            detector_names.push_back(pDet->Attribute("name"));
            detector_bool.push_back(strtod(pDet->Attribute("use"),&end));
            pDet = pDet->NextSiblingElement("detector");
            if(is_verbose)	std::cout<<"SBNconfig::SBnconfig\t|| loading detector: "<<detector_names.back()<<" with use_bool "<<detector_bool.back()<<std::endl;
        }
    }

    //How many channels do we want! At the moment each detector must have all channels
    int nchan = 0;
    if(!pChan){
        std::cout<<otag<<"ERROR: Need at least 1 channel defined in xml./n";
    }else{


        while(pChan){
            // Read in how many bins this channel uses
            channel_names.push_back(pChan->Attribute("name"));
            channel_units.push_back(pChan->Attribute("unit"));
            channel_bool.push_back(strtod(pChan->Attribute("use"),&end));
            num_bins.push_back(strtod(pChan->Attribute("numbins"), &end));

            if(is_verbose)	std::cout<<otag<<"Loading Channel : "<<channel_names.back()<<" with use_bool: "<<channel_bool.back()<<std::endl;

            // What are the bin edges and bin widths (bin widths just calculated from edges now)
            TiXmlElement *pBin = pChan->FirstChildElement("bins");
            std::stringstream iss(pBin->Attribute("edges"));

            double number;
            std::vector<double> binedge;
            std::vector<double> binwidth;
            while ( iss >> number ) binedge.push_back( number );

            for(int b = 0; b<binedge.size()-1; b++){
                binwidth.push_back(fabs(binedge.at(b)-binedge.at(b+1)));
            }
            if(binedge.size() != num_bins.back()+1){
                std::cout<<otag<<"ERROR: num_bins: "<<num_bins.back()<<" but we have "<<binedge.size()<<" binedges! should be num+1"<<std::endl;
                exit(EXIT_FAILURE);
            }

            bin_edges.push_back(binedge);
            bin_widths.push_back(binwidth);

            // Now loop over all this channels subchanels. Not the names must be UNIQUE!!
            TiXmlElement *pSubChan;

            pSubChan = pChan->FirstChildElement("subchannel");
            int nsubchan=0;
            while(pSubChan){
                //std::cout<<"Subchannel: "<<pSubChan->Attribute("name")<<" use: "<<pSubChan->Attribute("use")<<" osc: "<<pSubChan->Attribute("osc")<<std::endl;
                subchannel_names[nchan].push_back(pSubChan->Attribute("name"));
                subchannel_bool[nchan].push_back(strtod(pSubChan->Attribute("use"),&end));
                //0 means dont oscillate, 11 means electron disapearance, -11 means antielectron dis..etc..
                if(pSubChan->Attribute("osc"))
                {
                    has_oscillation_patterns = true;
                }

                subchannel_osc_patterns.at(nchan).push_back(strtod(pSubChan->Attribute("osc"), &end));

                if(is_verbose)	std::cout<<otag<<"--> Subchannel: "<<subchannel_names.at(nchan).back()<<" with use_bool "<<subchannel_bool.at(nchan).back()<<" and osc_pattern "<<subchannel_osc_patterns.at(nchan).back()<<std::endl;

                nsubchan++;
                pSubChan = pSubChan->NextSiblingElement("subchannel");
            }
            num_subchannels.push_back(nsubchan);

            nchan++;
            pChan = pChan->NextSiblingElement("channel");
        }
    }

    // if wea re creating a covariance matrix using a ntuple and weights, here is the info
    if(pMC){
        while(pMC)
        {
            //pot.push_back(strtof(pMC->Attribute("pot"),&end));
            //pot_scaling.push_back(strtof(pMC->Attribute("potscale"),&end));
            multisim_name.push_back(pMC->Attribute("treename"));
            multisim_file.push_back(pMC->Attribute("filename"));
            const char* maxevents = pMC->Attribute("maxevents");
            if(maxevents==NULL){
                multisim_maxevents.push_back(1e7);
            }else{
                multisim_maxevents.push_back(strtod(maxevents,&end) );
            }
            const char* scale = pMC->Attribute("scale");
            if(scale==NULL){
                multisim_scale.push_back(1.0);
            }else{
                multisim_scale.push_back(strtod(scale,&end) );
            }

            /*
            //Currently take all parameter variations in at once. Depreciated code. 
            TiXmlElement *pParams = pMC->FirstChildElement("parameters");
            std::stringstream sss(pParams->Attribute("names"));

            std::vector<std::string> vstring;
            std::string nam;
            while ( sss >> nam) vstring.push_back( nam );

            parameter_names.push_back(vstring);
            */

            TiXmlElement *pFriend;
            pFriend = pMC->FirstChildElement("friend");
            if(pFriend){

                if( multisim_file_friend_treename_map.count(multisim_file.back())>0){
                    (multisim_file_friend_treename_map[multisim_file.back()]).push_back( pFriend->Attribute("treename") );
                    (multisim_file_friend_map[multisim_file.back()]).push_back(pFriend->Attribute("filename"));
                }else{
                    std::vector<std::string> temp_treename = {pFriend->Attribute("treename")};
                    std::vector<std::string> temp_filename = {pFriend->Attribute("filename")};

                    multisim_file_friend_treename_map[multisim_file.back()] = temp_treename;
                    multisim_file_friend_map[multisim_file.back()] = temp_filename;
                }
            }


            TiXmlElement *pBranch;
            pBranch = pMC->FirstChildElement("branch");


            std::vector<BranchVariable*> TEMP_branch_variables;
            while(pBranch){

                std::string bnam = pBranch->Attribute("name");
                std::string btype = pBranch->Attribute("type");
                std::string bhist = pBranch->Attribute("associated_subchannel");

                //if(btype == "int"){
                //	std::cout<<"NO INT ALLOWED "<<bnam<<std::endl;
                //	exit(EXIT_FAILURE);
                //TEMP_branch_variables.push_back( new BranchVariable_i(bnam,btype, bhist ) );
                //}else
                if(btype == "double"){
                    if(is_verbose)std::cout<<otag<<"Setting double variable "<<bnam<<" @ "<<bhist<<std::endl;
                    TEMP_branch_variables.push_back( new BranchVariable_d(bnam, btype, bhist ) );
                    //}else if(btype == "float"){
                    //	std::cout<<otag<<"Setting float variable "<<bnam<<" @ "<<bhist<<std::endl;
                    //	TEMP_branch_variables.push_back( new BranchVariable_f(bnam, btype, bhist ) );
            }else{
                std::cout<<otag<<"ERROR: currently only double, allowed for reco variables (sorry!)\n";
                exit(EXIT_FAILURE);
            }

            std::string oscillate = "false";
            if(pBranch->Attribute("oscillate")!=NULL){
                oscillate =pBranch->Attribute("oscillate");
            }	

            if(oscillate == "false"){
                if(is_verbose)std::cout<<otag<<"Oscillations are Off. oscillate="<<oscillate<<std::endl;
                TEMP_branch_variables.back()->SetOscillate(false);
            }else if(oscillate == "true"){
                if(is_verbose)std::cout<<otag<<"Setting Oscillate! "<<oscillate<<std::endl;
                TEMP_branch_variables.back()->SetOscillate(true);
                TEMP_branch_variables.back()->true_param_name = pBranch->Attribute("true_param_name");
                TEMP_branch_variables.back()->true_L_name = pBranch->Attribute("true_L_name");
                std::cout<<otag<<"Set Oscillate! "<<pBranch->Attribute("true_param_name")<<" "<<pBranch->Attribute("true_L_name")<<std::endl;
            }else{
                if(is_verbose)std::cout<<otag<<"Do Not Oscillate "<<oscillate<<std::endl;
                TEMP_branch_variables.back()->SetOscillate(false);
            }


            pBranch = pBranch->NextSiblingElement("branch");
            }
            branch_variables.push_back(TEMP_branch_variables);
            //next file
            pMC=pMC->NextSiblingElement("MultisimFile");
        }
    }

    //Where is the "data" folder that keeps pre-computed spectra and rootfiles
    //Will eventuall just configure this in CMake
    while(pData){
        data_path = pData->Attribute("path");
        pData = pData->NextSiblingElement("data");
        if(is_verbose)std::cout<<otag<<"data path loaded as: "<<data_path<<std::endl;
    }


    // Where is the covariance matrix you want to use, and whats its name in the root file.
    while(pCov){
        correlation_matrix_rootfile = data_path + pCov->Attribute("file");
        correlation_matrix_name = pCov->Attribute("name");
        pCov = pCov->NextSiblingElement("covariance");
    }


    if(is_verbose) std::cout<<otag<<"|| Calculating used things"<<std::endl;

    // so num_channels here is number of TOTAL channels in xml.
    num_channels = channel_names.size();
    num_modes = mode_names.size();
    num_detectors  = detector_names.size();
    //Calculate bin_widths from bin_edges

    num_subchannels_xml = num_subchannels;
    num_channels_xml = num_channels;
    num_modes_xml = num_modes;
    num_detectors_xml = num_detectors;





    // here we run through every combination, and make note when (and where binwise) all the subchannels that are turned on are.
    std::string tempn;
    int indexcount = 0;

    for(int im = 0; im < num_modes; im++){

        for(int id =0; id < num_detectors; id++){

            for(int ic = 0; ic < num_channels; ic++){

                for(int sc = 0; sc < num_subchannels.at(ic); sc++){

                    tempn = mode_names[im] +"_" +detector_names[id]+"_"+channel_names[ic]+"_"+subchannel_names[ic][sc];
                    if(is_verbose)std::cout<<otag<<""<<tempn<<" "<<im<<" "<<id<<" "<<ic<<" "<<sc<<std::endl;

                    // This is where you choose NOT to use some fields
                    if(mode_bool[im] && detector_bool[id] && channel_bool[ic] && subchannel_bool[ic][sc]){

                        fullnames.push_back(tempn);
                        for(int k = indexcount; k < indexcount+num_bins.at(ic); k++){
                            used_bins.push_back(k);
                            //std::cout<<"USED: "<<k<<std::endl;

                        }
                    }
                    std::vector<int> tvec = {indexcount, indexcount+num_bins.at(ic)-1};

                    map_tag_to_covariance_index[tempn] = 	tvec;
                    indexcount = indexcount + num_bins.at(ic);


                }
            }
        }
    }





    //For here on down everything is derivable, above is just until I actually Get config working.
    if(is_verbose) std::cout<<otag<<"Starting on MC file parameters"<<std::endl;

    num_modes=0;
    for(int i=0;i<mode_bool.size(); i++){	if(mode_bool.at(i)){num_modes++; mode_used.push_back(i);}	}

    num_detectors = 0;
    for(int i=0; i<detector_bool.size(); i++){ if(detector_bool.at(i)) {num_detectors++; detector_used.push_back(i);}	}

    for(int i=0; i< num_channels; i++){
        num_subchannels.at(i) = 0;
        for(bool j: subchannel_bool[i]){ if(j) num_subchannels[i]++;}
    }
    //This needs to be above num_channel recalculation;

    num_channels = 0;
    for(int i=0; i< channel_bool.size(); i++){
        if( channel_bool.at(i)){
            num_channels++;
            channel_used.push_back(i);
        }
    }



    if(is_verbose) std::cout<<otag<<"Calculating Total Bins"<<std::endl;
    this->CalcTotalBins();


    if(is_verbose){
        std::cout<<otag<<"Checking number of XX"<<std::endl;
        std::cout<<otag<<"--> num_modes: "<<num_modes<<" out of "<<num_modes_xml<<std::endl;
        std::cout<<otag<<"--> num_detectors: "<<num_detectors<<" out of "<<num_detectors<<std::endl;
        std::cout<<otag<<"--> num_channels: "<<num_channels<<" out of "<<num_channels<<std::endl;
        for(auto i: channel_used){
            std::cout<<otag<<"----> num_subchannels: "<<num_subchannels.at(i)<<" out of "<<num_subchannels_xml.at(i)<<std::endl;
            std::cout<<otag<<"----> num_bins: "<<num_bins.at(i)<<std::endl;
        }

        std::cout<<otag<<"--> num_bins_detector_block: "<<num_bins_detector_block<<std::endl;
        std::cout<<otag<<"--> num_bins_detector_block_compressed: "<<num_bins_detector_block_compressed<<std::endl;
        std::cout<<otag<<"--> num_bins_mode_block: "<<num_bins_mode_block<<std::endl;
        std::cout<<otag<<"--> num_bins_mode_block_compressed: "<<num_bins_mode_block_compressed<<std::endl;


        std::cout<<otag<<"--> num_bins_total: "<<num_bins_total<<std::endl;
        std::cout<<otag<<"--> num_bins_total_compressed: "<<num_bins_total_compressed<<std::endl;


    }

    //Now delete all info corresponding to NON-USED channels & subchannels.


    //num_subchannels, num_bins, bin_edges, bin_widths, channel_names, subchannel_names;
    auto temp_num_subchannels = num_subchannels;
    auto temp_num_bins= num_bins;
    auto temp_subchannel_names = subchannel_names;
    auto temp_channel_names = channel_names;
    auto temp_bin_edges = bin_edges;
    auto temp_bin_widths = bin_widths;
    auto temp_detector_names = detector_names;
    auto temp_mode_names = mode_names;
    auto temp_mode_bool = mode_bool;
    auto temp_channel_bool= channel_bool;
    auto temp_detector_bool = detector_bool;
    auto temp_subchannel_bool = subchannel_bool;
    auto temp_subchannel_osc_patterns = subchannel_osc_patterns;

    num_subchannels.clear();
    num_bins.clear();
    subchannel_names.clear();
    channel_names.clear();
    bin_edges.clear();
    bin_widths.clear();
    detector_names.clear();
    mode_names.clear();

    mode_bool.clear();
    channel_bool.clear();
    subchannel_bool.clear();
    detector_bool.clear();

    subchannel_osc_patterns.clear();


    for(int c: channel_used){
        //if(is_verbose){std::cout<<otag<<"Adding channel: "<<c<<std::endl;}
        num_subchannels.push_back( temp_num_subchannels.at(c));
        num_bins.push_back( temp_num_bins.at(c));
        subchannel_names.push_back( temp_subchannel_names.at(c));
        channel_names.push_back( temp_channel_names.at(c));
        bin_edges.push_back( temp_bin_edges.at(c));
        bin_widths.push_back( temp_bin_widths.at(c));

        channel_bool.push_back(temp_channel_bool.at(c));
        subchannel_bool.push_back(temp_subchannel_bool.at(c));
        subchannel_osc_patterns.push_back(temp_subchannel_osc_patterns.at(c));
    }
    for(int d: detector_used){
        detector_names.push_back(temp_detector_names.at(d));
        detector_bool.push_back(temp_detector_bool.at(d));
        //if(is_verbose) std::cout<<otag<<"Using Detector: "<<detector_names.back()<<std::endl;
    }

    for(int m: mode_used){
        mode_names.push_back(temp_mode_names.at(m));
        mode_bool.push_back(temp_mode_bool.at(m));
    }



    if(is_verbose) {

        //	std::cout<<"--> num_channels: "<<num_channels<<" channel_bool.size(): "<<channel_bool.size()<<" channel_names.size(): "<<channel_names.size()<<std::endl;
        //	std::cout<<"--> num_modes: "<<num_modes<<" mode_bool.size(): "<<mode_bool.size()<<" mode_names.size(): "<<mode_names.size()<<std::endl;
        //	std::cout<<"--> num_detectors: "<<num_detectors<<" detector_bool.size(): "<<detector_bool.size()<<" detector_names.size(): "<<detector_names.size()<<std::endl;
        for(int i=0; i< num_channels; i++){
            //		std::cout<<"--> num_subchannels: "<<num_subchannels.at(i)<<" subchannel_bool.size(): "<<subchannel_bool.at(i).size()<<" subchannel_names.at(i).size(): "<<subchannel_names.at(i).size()<<std::endl;
        }
    }



    a_num_bins = num_bins.data();
    a_num_subchannels = num_subchannels.data();



    if(is_verbose){std::cout<<otag<<"Done!"<<std::endl;}




}//end constructor

SBNconfig::SBNconfig(std::string whichxml): SBNconfig(whichxml, true) {}

//Constructor to build a SBNspec from scratch, not reccomended often
SBNconfig::SBNconfig(std::vector<std::string> modein, std::vector<std::string> detin, std::vector<std::string> chanin, std::vector<std::vector<std::string>> subchanin, std::vector<std::vector<double>> binin){

    otag = "SBNconfig::SBNconfig\t|| ";
    is_verbose = true;

    num_detectors = detin.size();
    num_channels = chanin.size();
    num_modes = modein.size();

    if(subchanin.size() != chanin.size()){
        std::cout<<otag<<"ERROR SUBCHAN.size() != chanin.size()"<<std::endl;
        exit(EXIT_FAILURE);
    }

    for(auto sb: subchanin){
        num_subchannels.push_back( sb.size());
    }

    for(auto bn: binin){
        num_bins.push_back(bn.size()-1);
    }

    this->CalcTotalBins();

    xmlname = "NULL";


    //the xml names are the way we track which channels and subchannels we want to use later
    mode_names = modein;
    detector_names = detin;
    channel_names = chanin;
    subchannel_names = subchanin;

    for(auto c: chanin){
        channel_bool.push_back(true);
    }
    for(auto m: modein){
        mode_bool.push_back(true);
    }
    for(auto d: detin){
        detector_bool.push_back(true);
    }
    for(auto c: chanin){
        std::vector<bool> tml;
        for(int i=0; i< num_subchannels.size(); i++){
            tml.push_back(true);
        }
        subchannel_bool.push_back(tml);
    }

    //self explanatory
    bin_edges = binin;

    for(auto binedge: bin_edges){
        std::vector<double> binwidth;
        for(int b = 0; b<binedge.size()-1; b++){
            binwidth.push_back(fabs(binedge.at(b)-binedge.at(b+1)));
        }
        bin_widths.push_back(binwidth);
    }

    // The order is IMPORTANT its the same as defined in xml
    for(auto m: mode_names){
        for(auto d: detector_names){
            for(int c = 0; c< num_channels; c++){
                for(auto sb: subchannel_names.at(c)){
                    std::string tmp = m +"_"+ d+"_"+channel_names.at(c)+"_"+sb;
                    fullnames.push_back(tmp);
                }
            }
        }
    }


    for(int i=0; i< num_bins_total; i++){
        used_bins.push_back(i);
    }


};


int SBNconfig::CalcTotalBins(){

    // These variables are important
    // They show how big each mode block and decector block are, for any given number of channels/subchannels
    // both before and after compression!

    //needs to be calculated AFTER usage bool removal above

    num_bins_detector_block = 0;
    num_bins_detector_block_compressed = 0;

    for(auto i: channel_used){
        num_bins_detector_block += num_subchannels[i]*num_bins[i];
        num_bins_detector_block_compressed += num_bins[i];
    }

    num_bins_mode_block = num_bins_detector_block*num_detectors;
    num_bins_mode_block_compressed = num_bins_detector_block_compressed*num_detectors;

    num_bins_total = num_modes*num_bins_mode_block;
    num_bins_total_compressed = num_modes*num_bins_mode_block_compressed;

    return 0;
}
