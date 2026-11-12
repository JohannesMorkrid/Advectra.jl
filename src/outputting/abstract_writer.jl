# ------------------------------------------------------------------------------------------
#                                 Abstract Writer Interface
# ------------------------------------------------------------------------------------------

"""
    AbstractWriter

  Abstract supertype for all output writers. A writer is responsible for storing diagnostic
  samples produced during a simulation. Concrete subtypes decide *where* and *how* data is
  persisted (HDF5 file, in-memory arrays, NetCDF, Zarr, …).

  # Interface

  ## Required — every concrete subtype must implement:

  - `initialize!(writer, prob, diagnostics, initial_samples, t0; kwargs...)`:
      Called once before the time loop by `Output`. Should build the file structure and
      write the initial sample for each diagnostic. Currently recognised kwargs passed
      by `Output`: `simulation_name`, `resume`. Writers extract only what they need and
      absorb the rest via `kwargs...`.

  - `write_to!(writer, diagnostic, sample, step, time)`:
      Called every time a diagnostic is sampled. Should append `sample` to the writer's
      storage. `step` serves as a commit marker for crash safety. Writers determine their
      own write position — no `idx` is passed. Implementations may check
      `isinitialized(writer)` and throw if called before `initialize!`.

  ## Optional — default no-ops provided:

  - `flush!(writer, step)`:
      Called after every step. May be a no-op for writers that do not buffer (e.g.
      `LocalWriter`).

  - `Base.close(writer)`:
      Called at the end of a simulation. Should release any open file handles or resources.

  - `Base.isopen(writer)`:
      Returns whether the writer is open and ready to accept writes. Guarding against writes
       after `close` has been called.

  - `store_checkpoint!(writer, cache, step, time)`:
      Store a checkpoint of the integrator `cache` at `step`/`time`. Defaults to no-op.
      
  ## Not given a default — must be implemented when `supports_checkpoint` returns `true`:

  - `restore_checkpoint(writer, prob, scheme) -> (cache, step, t)`:
      Restore the integrator cache, step and time from a previously stored checkpoint.
      Only called on the first writer for which `has_checkpoint` returns `true`. Has no
      default no-op because silently returning nothing would leave the simulation in an
      uninitialised state.

  ## Trait predicates:

  - `isinitialized(writer)`:
      Returns `true` once `initialize!` has been called. Used to guard `write_to!` calls.

  - `supports_checkpoint(writer)`:
      Returns `true` for writers that implement `store_checkpoint!` and
      `restore_checkpoint`. Used by `Output` to find the restore source.

  - `has_checkpoint(writer)`:
      Returns `true` when a checkpoint has actually been stored and can be restored.
      Distinct from `supports_checkpoint` — a writer may support checkpointing but not
      yet have one stored (e.g. a fresh run).
"""
abstract type AbstractWriter end

# ------------------------------- Default no-op implementations ---------------------------

"""Default: writers that do not buffer have nothing to flush."""
flush!(::AbstractWriter, step) = nothing

"""Default: writers with no file handle have nothing to close."""
Base.close(::AbstractWriter) = nothing

"""Default: writers are considered open unless overridden."""
Base.isopen(::AbstractWriter) = true

"""Default: writers that do not support checkpointing silently skip."""
store_checkpoint!(::AbstractWriter, cache, step, time) = nothing

"""Default: writers do not support checkpointing unless overridden."""
supports_checkpoint(::AbstractWriter) = false

"""Default: writers have no stored checkpoint unless overridden."""
has_checkpoint(::AbstractWriter) = false

"""Default: writers are considered uninitialized until overridden."""
isinitialized(::AbstractWriter) = false