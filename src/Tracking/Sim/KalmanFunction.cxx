#include "Tracking/Sim/CKFProcessor.h"

using namespace Acts::UnitLiterals;

namespace tracking {
namespace sim {

void CKFProcessor::kalmanRefit(const Acts::MultiTrajectory &mj,
                          std::size_t trackTip,
                          const LdmxMeasurementCalibrator &calibrator,
                          const Acts::PropagatorOptions<ActionList, AbortList> &prop_options,
                          const Acts::Surface *target_surface,
                          const Acts::BoundTrackParameters &start_parameters,
                          Acts::Logging::Level log_level) {
  //   std::cout<<"Preparing theKF refit"<<std::endl; 
  std::vector<std::reference_wrapper<const ActsExamples::IndexSourceLink>> fit_trackSourceLinks;
  mj.visitBackwards(trackTip, [&](const auto& state) {
    if( state.hasUncalibrated() ) {
      const auto& sourceLink =
        static_cast<const ActsExamples::IndexSourceLink&>(state.uncalibrated());
      auto typeFlags = state.typeFlags();
      if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
        fit_trackSourceLinks.push_back(std::cref(sourceLink));
      }
    }
  });

  const auto kfLogger = Acts::getDefaultLogger("KalmanFitter", log_level);
  
  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;
  Acts::KalmanFitterExtensions kfitter_extensions;
  kfitter_extensions.calibrator.connect<&LdmxMeasurementCalibrator::calibrate_1d>(&calibrator);
  kfitter_extensions.updater.connect<&Acts::GainMatrixUpdater::operator()>(&kfUpdater);
  kfitter_extensions.smoother.connect<&Acts::GainMatrixSmoother::operator()>(&kfSmoother);

  //rFiltering is true, so it should run in reversed direction.
  Acts::KalmanFitterOptions kfitter_options{gctx_,bctx_,cctx_,
                                kfitter_extensions,Acts::LoggerWrapper{*kfLogger},
                                prop_options,target_surface,
                                true, true, true}; //mScattering, exoLoss, rFiltering
  


  
  // create the Kalman Fitter
  if (debug_) {
    std::cout<<"Make the KalmanFilter fitter object"<<std::endl;
    std::cout<<"Refit"<<std::endl;
    std::cout<<"Starting from "<<std::endl;
    std::cout<<start_parameters.position(gctx_)<<std::endl;
    std::cout<<"With momenutm"<<std::endl;
    std::cout<<start_parameters.momentum()<<std::endl;
    std::cout<<"rFiltering =" << std::boolalpha << kfitter_options.reversedFiltering <<std::endl;
  }
  
  auto kf_refit_result = kf_->fit(fit_trackSourceLinks.begin(),fit_trackSourceLinks.end(),
                                start_parameters, kfitter_options);
  
  if (not kf_refit_result.ok()) {
    std::cout<< "KF Refit failed in event " << nevents_ << std::endl;
    return;
  } 
  
  auto kf_refit_value = kf_refit_result.value();
  auto kf_params = kf_refit_value.fittedParameters;
  
  if (not kf_params) {
      std::cout << "KF result does not contain fitted parameters\n";
      return;
  }
  
  if ( kf_params->absoluteMomentum() > 4.5_GeV ) {
    std::cout << "KF refit momentum is " << kf_params->absoluteMomentum() << " GeV in event " << nevents_ << "\n";
  }
  
  h_p_refit_     ->Fill((*kf_params).absoluteMomentum());
  h_d0_refit_    ->Fill((*kf_params).get<Acts::BoundIndices::eBoundLoc0>());
  h_z0_refit_    ->Fill((*kf_params).get<Acts::BoundIndices::eBoundLoc1>());
  h_phi_refit_   ->Fill((*kf_params).get<Acts::BoundIndices::eBoundPhi>());
  h_theta_refit_ ->Fill((*kf_params).get<Acts::BoundIndices::eBoundTheta>());
}

}
}
