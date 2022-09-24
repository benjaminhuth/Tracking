#include "Tracking/Sim/CKFProcessor.h"

using namespace Acts::UnitLiterals;

namespace tracking {
namespace sim {

void CKFProcessor::gsfRefit(const Acts::MultiTrajectory &mj,
                            std::size_t trackTip,
                            const LdmxMeasurementCalibrator &calibrator,
                            const Acts::PropagatorOptions<ActionList, AbortList> &prop_options,
                            const Acts::Surface *target_surface,
                            const Acts::BoundTrackParameters &start_parameters,
                            Acts::Logging::Level log_level) {
  const auto gsfLogger = Acts::getDefaultLogger("GSF", log_level);
  std::vector<std::reference_wrapper<const ActsExamples::IndexSourceLink>> fit_trackSourceLinks;
  mj.visitBackwards(trackTip, [&](const auto& state) {
    auto typeFlags = state.typeFlags();
    if (typeFlags.test(Acts::TrackStateFlag::MeasurementFlag)) {
      const auto& sourceLink =
          static_cast<const ActsExamples::IndexSourceLink&>(state.uncalibrated());
      fit_trackSourceLinks.push_back(std::cref(sourceLink));
    }
  });

  //Same extensions of the KF
  Acts::GainMatrixUpdater kfUpdater;
  Acts::GsfExtensions gsf_extensions;
  gsf_extensions.calibrator.connect<&LdmxMeasurementCalibrator::calibrate_1d>(&calibrator);
  gsf_extensions.updater.connect<&Acts::GainMatrixUpdater::operator()>(&kfUpdater);

  Acts::GsfOptions gsf_options{gctx_,
    bctx_,
    cctx_,
    gsf_extensions,
    Acts::LoggerWrapper{*gsfLogger},
    prop_options,
    target_surface};

  gsf_options.abortOnError = false;
  gsf_options.maxComponents = gsfComponents_;
  gsf_options.disableAllMaterialHandling = false;
  
  auto gsf_refit_result = gsf_->fit(fit_trackSourceLinks.begin(),
                                    fit_trackSourceLinks.end(),
                                    start_parameters,
                                    gsf_options);
  
  if (not gsf_refit_result.ok()) {
    std::cout << "GSF Refit failed for event " << nevents_ << std::endl;
    return;
  }
  
  auto gsf_refit_value = gsf_refit_result.value();
  auto gsf_params = gsf_refit_value.fittedParameters;
  
  if ( not gsf_params ) {
    std::cout << "GSF result does not contain fitted parameters\n";
    return;
  }
  
  if ( gsf_params->absoluteMomentum() > 4.5_GeV ) {
    std::cout << "GSF refit momentum is " << gsf_params->absoluteMomentum() << " GeV in event " << nevents_ << "\n";
  }
  
  h_p_gsf_refit_     ->Fill(gsf_params->absoluteMomentum());
  h_d0_gsf_refit_    ->Fill(gsf_params->get<Acts::BoundIndices::eBoundLoc0>());
  h_z0_gsf_refit_    ->Fill(gsf_params->get<Acts::BoundIndices::eBoundLoc1>());
  h_phi_gsf_refit_   ->Fill(gsf_params->get<Acts::BoundIndices::eBoundPhi>());
  h_theta_gsf_refit_ ->Fill(gsf_params->get<Acts::BoundIndices::eBoundTheta>());
}

}
}
