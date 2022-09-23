{
/*

gROOT->ProcessLine(".L Hypo2ChPionsPi02Photons.cpp++");
gROOT->ProcessLine(".L Hypo2ChPions2PhotonsLostPi0.cpp++");
gROOT->ProcessLine(".L Hypo2ChPions2pi0.cpp++");
*/
gROOT->ProcessLine(".L Hypo2ChPionsPi02Photons_cpp.so");
gROOT->ProcessLine(".L Hypo2ChPions2PhotonsLostPi0_cpp.so");
gROOT->ProcessLine(".L Hypo2ChPions2pi0_cpp.so");


gROOT->ProcessLine(".L TrPh.C++");
gROOT->ProcessLine("callLoop()");

}
