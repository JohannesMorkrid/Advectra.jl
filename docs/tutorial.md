# Tutorial
## File structure overview
The *src* folder holds the code and modules containing the structs (objects) and methods used in this project. Of note are the ``Domain``, ``SpectralODEProblem`` and ``Output`` structs. These structs are used to initialize the simulations and will be explained below. The most important method is the *spectral_solve* method used to evolve the solution in time and sample the data. The *tests* folder houses the different simulations ran, the *code testing* used for code validation and some code used for development in the *development* folder.

The folder of interest is the *tests/Sheath-interchange/simulation batch* where the [gyro-bohm=\<sigma value\>.jl](<../tests/Sheath-interchange/simulation batch/gyro-bohm=1e-1.jl>) files are the once used to get the data analyzed in section 5 of the thesis. In addition the *tests/Garcia 2006 Pop* and *tests/Kube 2011 Pop* folders were used for blob simulations.

The *docs* folder is for documentation, which houses this tutorial, but was downprioritized during the master project writing.

## Importing the modules
As this code is not ''Released'' the following code needs to be included to import the code into a ``.jl`` file:
```
using Advectra
```
**NOTE:** This assumes the user is in the *Advectra* folder.

## Running a simulation
To run a simulation the following is needed:
* A ``SpectralODEProblem`` which contains the linear operator ``L``, the non-linear operator ``N``, the ``Domain`` which contains info about the domain used by the operators, an initial condition (abreviated ic) in the form of an ``AbstractArray`` often a $(N_y,N_x,2)$ ``Matrix``, a timespan specifying start $t_0$ and stop times $t_{\textrm{end}}$ (tuple with 2 entries), a parameters ``Dict`` (of the same data type) and a timestep $dt$.
* A scheme, of which $\textrm{MSS}p$ up to order $p=3$ is the only supported. Note the scheme has to be initialized in the function call (for instance ``MSS3()``) to get the right ``DataType``.
* An optional ``Output`` struct, which containts information about the diagnostics, sampling rates and naming conventions. If no ``Output`` is passed, a file with a random filename is generated storing only the output data.

In addition the code supports the option to ''resume'' the code which is toggled using the ``resume`` keyword argument. To interupt the program use ``Ctrl+C`` (currently buggy if not ran in interactive mode).

## SpectralOperators
The code is written in such a way that the spatial derivatives uses the following Template ``operatorName(fields..., domain)``. These spatial derivatives are refered to as SpectralOperators in the code and uses ``SpectralOperatorCache`` stored in the ``Domain`` to allocate coefficients used by the Linear operators, and cached ``AbstractArray``'s used mostly by the Non-linear operators. These operators are *in-place* and returns the altered array.

The linear operator method ``L`` is used to calculate the scaling coefficient $c$ in the $\textrm{MSS}p$ scheme, with the output assumed constant, while the non-linear operator ``N`` is used in the explicit part of the scheme.

The non-linear operator ``N`` looks like:
```
function N(u, d, p, t)
    n = @view u[:, :, 1]
    Ω = @view u[:, :, 2]
    ϕ = solve_phi(Ω, d)

    dn = -poisson_bracket(ϕ, n, d)
    dn .-= (1 - p["g"]) * diff_y(ϕ, d)
    dn .-= p["g"] * diff_y(n, d)
    dn .+= p["sigma"] * ϕ

    dΩ = -poisson_bracket(ϕ, Ω, d)
    dΩ .-= p["g"] * diff_y(n, d)
    dΩ .+= p["sigma"] * ϕ
    return cat(dn, dΩ, dims=3)
end
```
where it is worth noting the use of the ``@view`` macro to access a slice of the solution $u$ (in spectral space so should really be $\hat{u}$), the ``solve_phi`` used to obtain the electrostatic potential from the vorticity field, the ``.`` infront of the operators ensure the operators happens elementwise and also speeds up the calculation and ``cat(dn, dΩ, dims=3)`` is used to combine the two changes to the field.

The operators supported are ``diff_x``, ``diff_xx``, ``diff_y`` and ``diff_yy`` for spatial derivatives in the $x$ and $y$ direction respectively. The ``laplacian`` operator with aliases ``diffusion`` and ``Δ`` with the support of ``hyper_diffusion``. ``solve_phi`` which computes the electrostatic potential from the vorticity. The non-linear operators supported are ``quadratic_term`` which computes the quadratic term using the pseudo-spectral method, ``poisson_bracket`` which utilizes the ``quadratic_term`` method, and also ``reciprocal``, ``spectral_exp``, ``spectral_expm1`` and ``spectral_log`` which all uses ``spectral_function`` to evaluate the function in physical space and get the modes back.

## Domain
The ``Domain`` contains the information about the spatial discretization, the ``SpectralOperatorCache`` (explained above) and the transform plans (``FFTPlans`` in the current implementation). The easiest way to initialize a domain is ``Domain(N)`` which assumes a centered unit square domain. One may also use ``Domain(N,L)`` for a square domain with a different size. Both initilizations are optimized for run-time, for full flexibility use
```
Domain(Nx, Ny, Lx, Ly; real_transform=true, dealiased=true, use_cuda=true, x0=-Lx / 2, y0=-Ly / 2)
```
where ``x0`` and ``y0`` specify the down-left corner position (origin). The keyword arguments allows for the user to toggle the use of de-aliasing ``anti-aliasing``, realFFT ``real_transform`` and the use of *CUDA* ``use_cuda``,currently only supported in the *cuda* branch. Based on the toggles the program allocates fitting FFTPlans (the ``transform`` field of the ``Domain``) and ``SpectralOperatorCache``. The ``Domain`` discretization is stored in the ``x`` and ``y`` fields with the respective wavenumbers stored in ``kx`` and ``ky``.

## Initial condition
The initial condition needs to be compatible with the expected return value of the linear and non-linear operators. To achieve this use the ``initial_condition(fun, domain)`` function for a specified function ``fun`` on the form ``f(x,y;kwargs...)`` or use the badly named ``initial_condition_linear_stability(domain, amplitude)`` for random uncorrelated modes with a cross-phase of $\pi/2$ and amplitude determined by ``amplitude``.

## SpectralODEProblem
The ``SpectralODEProblem`` is initiated using 
```
SpectralODEProblem(L, N, u0, domain, tspan; p=Dict(), dt=0.01, remove_modes=remove_nothing, kwargs...
```
where ``L`` and ``N`` are the aformentioned operators, ``u0`` is the initial condition, ``domain`` is of type ``Domain`` and ``tspan`` is the time span/interval $[t_0, t_{\textrm{end}}]$. A parameter dictionary ``p`` can be passed into the operators to allow for the use of parameter values. The ``SpectralODEProblem`` also contains the method for removing certain modes, as the initial condition ``u0`` in physical space is transformed into the spectral counterpart ``u0_hat``, on which the removal of mode is also imposed. Additional key word arguments ``kwargs`` may be stored in the ``SpectralOperatorCache`` to be passed to the ``Diagnostic``'s.

## Output and Diagnostics
The ``Output`` struct is responsible for sampling the ``Diagnostic``'s specified in the ``diagnostic`` 'list' and outputting it to a ``.h5`` file. To initilize the output use
```
Output(prob, N_samples, diagnostics, filename=basename(tempname())*"h5"; physical_transform=identity, simulation_name=:timestamp, store_hdf=true, store_locally=true, h5_kwargs...)
```
The ``N_samples`` is the number of raw solutions (data) outputed to the file. If ``filename`` is not specified a random string will be generated. A file can hold multiple ``simulation``'s (field of the ``Output`` struct), and the ``simulation_name`` can be user defined, the ``:timestamp`` on startup or a string formed by the ``:parameters``. A ``physical_transform`` can be applied to the solution in physical space before applying diagnostics or storing the fields, such as applying the exponential to the 'density field'. The ``Output`` struct also has the toggles to store the data locally (``store_locally``) and or store the data to a *HDF5* file (``store_hdf``).

Diagnostics is detailed in the master thesis. # TODO add to documentation!

## Additional functionalitites
The ``send_mail`` method allows the code to send a mail to the user once a simulation is finnished, but requires the user to configure the *Advectra/.env* file. Send a mail to get the app password for simulation.update@gmail.com.