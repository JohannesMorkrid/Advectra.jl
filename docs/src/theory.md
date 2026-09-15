## State
The state $S\in \mathcal{D}$ is a collection of fields $f_i$, where i=1, 2, ..., N_d$, expressed as
$$S(\bm{r},t) = [f_1(\bm{r},t), f_2(\bm{r},t), ..., f_d(\bm{r},t)].$$
It is a contiguous tensor/Array that can undergo standard time evolution scheme operations, such as elementwise operations, and which also can be unpacked into its individual fields using either a accesor pattern $S.f_i$ or using a the `unpack_state` method. It is also usefull to map a operation to each field using the `map_field[!]` operator.

## Field
A field $f_i\in\mathcal{D}$ is a view into a contiguous slice of the state $S$. The benefit of using a field over an Array is the fact that the domain $\mathcal{D}$ is attached and can easily be accessed using the `get_domain` method or through the `f_i.domain` access. Thus informing the spectral operators about the domain the field is defined on, the transform methods as well as the supported wave vectors.