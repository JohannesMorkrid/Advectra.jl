# ------------------------------------------------------------------------------------------
#                                    Abstract Writer Interface
# ------------------------------------------------------------------------------------------

"""
    AbstractWriter

  Abstract supertype for all output writers. A writer is responsible for storing diagnostic
  samples produced during a simulation. Concrete subtypes decide *where* and *how* data is
  persisted (HDF5 file, in-memory arrays, NetCDF, Zarr, …).

  # Interface — every concrete subtype must implement:

  - `initialize!(writer, prob, diagnostics, initial_samples, strides, t0)`:
      Called once before the time loop. Should allocate storage and write the initial sample.

  - `write!(writer, diagnostic, idx, sample, time)`:
      Called every time a diagnostic is sampled. Should persist `sample` at index `idx`.

  - `flush!(writer, step)`:
      Called after every step. May be a no-op for writers that do not buffer.

  - `Base.close(writer)`:
      Called at the end of a simulation. Should release any open file handles or resources.
"""
abstract type AbstractWriter end

# Default no-op implementations so writers only have to override what they need.
flush!(::AbstractWriter, step) = nothing
import Base.close
Base.close(::AbstractWriter) = nothing
