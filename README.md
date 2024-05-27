# ICB
Updated configuration for testing the iceberg module in NEMO

NEMO 4.2 offered an idealized configuration for testing iceberg dynamics and thermodynamics called ICB. Its description suggested that the dynamics and thermodynamics of the liquid ocean domain were not allowed to evolve in time. That was supposed to provide a suitable environment for testing iceberg code only.

However, at some point NEMO switched to modified leapfrog timestepping. As a result, the ocean dynamics and thermodynamics in the ICB test configuration were no longer fixed in time.

The updated ICB configuration provided here allows the NEMO user to hold the ocean state constant while testing iceberg code.
