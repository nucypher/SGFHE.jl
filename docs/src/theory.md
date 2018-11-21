# Theory

This section contains an extract from [S. Gao, "Efficient fully homomorphic encryption scheme"](https://eprint.iacr.org/2018/637) with the algorithms used in the package. For the extended discussion and proofs please refer to the original paper.

The algorithms described here have some minor changes as compared to the original paper, which preserve the mathematical idea behind them, but make implementation more straightforward. For example, all operations in this section are performed on unsigned integers, while the original paper has some expressions where signed and unsigned (modulo) integers are mixed and the conversion is performed implicitly. There are some minor changes in notation as well.


## Notation

``\mathbb{Z}_q = \mathbb{Z}/q\mathbb{Z}`` is the ring of integers modulo ``q``, and ``R_{n,q} = \mathbb{Z}[x]/(x^n+1, q)`` is a ring of polynomials modulo ``x^n+1`` (that is, negacyclic) with the coefficients from ``\mathbb{Z}_q``.

Vectors are denoted by bold-style symbols: ``\boldsymbol{a}`` and treated as row vectors. Their elements are accessed as ``a_i``. Matrices are denoted by capital letters: ``A``. Polynomials are denoted by writing their parameter explicitly: ``a(x)``.

A polynomial can be created out of a vector of coefficients with the operator ``\mathcal{P}[\boldsymbol{a}] = \sum_{i=0}^{n-1} a_i x^i``, where ``\boldsymbol{a} \in \mathbb{Z}_q^n``. By default, the result is belongs to ``R_{n,q}``, although it may be treated as belonging to some ``R_{m,Q}`` with ``m \ge n`` and ``Q \ge q``, in which case it will be stated explicitly.

Conversely, a vector of coefficients can be created out of a polynomial with the operator ``\mathcal{C}[a(x)] = \mathcal{C}[\mathcal{P}[\boldsymbol{a}]] = a_i``.

``x \div y \equiv \lfloor x/y \rfloor`` denotes floor division of unsigned integers.

``\mathcal{U}(a, b)`` denotes a uniform random integer from the range ``[a, b]``.


### Negative powers in polynomials

We make an occasional use of polynomials where some terms have negative powers (see [Bootstrap](@ref Bootstrap-theory) and [Packing LWEs](@ref)), that is linear combinations of expressions ``x^{-j} \,\mathrm{mod}\,(x^n+1)``. The modulus in this case is formally taken as
```math
x^{-j} \,\mathrm{mod}\,(x^n+1)
= x^{-j} - x^{-j} (x^n+1)
= -x^{n-j}.
```

In other words, the "left shift" operation (division by ``x^j``) works the same way as the "right shift" operation (multiplication by ``x^j``): if the term crosses the boundary of ``0`` or ``n-1``, respectively, its coefficient changes sign.


## Helper functions

This subsection describes some utility functions used in the FHE scheme.


### Pseudo-random expansion (``\mathrm{Expand}``)

*Paper:* Section 4.1.

*Implementation:* [`prng_expand()`](@ref SGFHE.prng_expand)

Pseudo-random expansion ``\mathrm{Expand}(\boldsymbol{u}, t)`` where ``\boldsymbol{u} \in \{0, 1\}^n`` produces a pseudo-random bit matrix ``\boldsymbol{V} \in \{0, 1\}^{n \times t}``, using ``\boldsymbol{u}`` as the seed. The result is converted to ``\boldsymbol{a} \in \mathbb{Z}^n`` by treating each row as the bit representation of an unsigned integer:
```math
a_i = \sum_{j=0}^t V_{ij} 2^j, \quad i = 0, \dots, n-1.
```

See the paper for the discussion about the properties the PRNG must have and possible implementations.


### Subvector extraction (``\mathrm{Extract}``)

*Paper:* Section 3.

*Implementation:* [`extract()`](@ref SGFHE.extract)

The function ``\mathrm{Extract}[\boldsymbol{a}, i, n]``, where ``\boldsymbol{a} \in \mathbb{Z}_q^m``, ``m \ge n``, and ``0 \le i < m`` is defined as
* if ``i < n - 1``,
  ```math
  \mathrm{Extract}[\boldsymbol{a}, i, n] = (a_i, a_{i-1}, \dots, a_0, -a_{m-1}, -a_{m-2}, \dots, -a_{m-1-i}),
  ```
* otherwise,
  ```math
  \mathrm{Extract}[\boldsymbol{a}, i, n] = (a_i, a_{i-1}, \dots, a_{i-n+1}).
  ```

The result is a vector in ``\mathbb{Z}_q^n``.


### Modulus reduction (``\mathrm{ModRed}``)

*Paper:* Lemma 2.3.

*Implementation:* [`reduce_modulus()`](@ref SGFHE.reduce_modulus)

The function ``\mathrm{ModRed}[a, q_1, q_2]`` where ``0 \le x < q_1`` and ``q_2 < q_1`` scales ``a \in \mathbb{Z}_{q_1}`` to a new modulus proportionally, with the result being in ``\mathbb{Z}_{q_2}``:
```math
\mathrm{ModRed}[a, q_1, q_2] = \left\lfloor \frac{a q_2}{q_1} \right\rceil \,\mathrm{mod}\,q_2.
```

``\mathrm{ModRed}`` can be applied to a vector, converting each of its elements. That is, if ``\boldsymbol{b} = \mathrm{ModRed}[\boldsymbol{a}, q_1, q_2]``, then the elements of ``\boldsymbol{b}`` are ``b_i = \mathrm{ModRed}[a_i, q_1, q_2]``.

Furthermore, ``\mathrm{ModRed}`` can be applied to the objects of belonging to ``\mathbb{Z}_q^n \times \mathbb{Z}_q`` (in particular, [LWE ciphers](@ref Creating-LWEs)), as
```math
\mathrm{ModRed}[(\boldsymbol{a}, b), q_1, q_2]
= (\mathrm{ModRed}[\boldsymbol{a}, q_1, q_2], \mathrm{ModRed}[b, q_1, q_2]).
```


### [Flattening (``\triangleleft``)](@id Flattening)


#### Integer flattening

*Paper:* Section 2.5.

*Implementation:* [`flatten()`](@ref SGFHE.flatten)

Given ``a \in \mathbb{Z}_q``, and positive integers ``B`` and ``\ell``, ``B^\ell \ge q``, ``a`` can be uniquely decomposed into a vector ``\boldsymbol{u} \in \mathbb{Z}_q^\ell`` such that
```math
a = \sum_{j=0}^{\ell-1} u_j B^j \,\mathrm{mod}\,q,
```
where ``0 \le u_j < B``. This can be done simply by successive division of ``a`` by ``B^j``, ``j = \ell-1 \dots 0``, storing the result as ``u_j`` and substituting the remainder for ``a`` in the next iteration.

Since in the above equation the right side is calculated modulo ``q``, we can shift the range of the elements of ``\boldsymbol{u}`` to any interval of length ``B``: ``s \le u_j < s + B``. This can be done by first adding the offset to ``a``:
```math
a^\prime = a + s \sum_{j=0}^{\ell-1} B^j \,\mathrm{mod}\,q.
```
Then the resulting ``a^\prime`` is flattened as described above, into ``\boldsymbol{u}^\prime``. Finally, the offset is subtracted from each element of the result:
```math
u_j = u_j^\prime - s \,\mathrm{mod}\,q.
```

In this FHE scheme the interval is chosen to be ``-B/2 < u_j \le B/2`` (and since we are dealing with non-negative numbers, comparisons undestood modulo ``q``, that is ``u_j \le B/2`` or ``u_j > q - B/2``), so in the algorithm above ``s = B/2-1`` (``-1`` to account for the left side of the interval being open instead of closed).

We denote this operation ``\boldsymbol{u} = a \triangleleft (B, \ell)``.


#### Random integer flattening

*Paper:* Lemma 2.4.

*Implementation:* [`flatten()`](@ref SGFHE.flatten)

The operation above can be randomized by allowing the result's elements to lie in a larger interval, in our case ``-2B < u_j \le 2B`` (with comparisons again understood modulo ``q``):

* Generate ``v_j = \mathcal{U}(-3B/2, 3B/2)``, ``j = 0, \dots, \ell-1``.
* Flatten deterministically ``\boldsymbol{u}^\prime = (a - \sum_{j=0}^{\ell-1} v_j B^j \,\mathrm{mod}\,q) \triangleleft (B, \ell)``.
* The result ``\boldsymbol{u}`` has elements ``u_j = u^\prime_j + v_j \,\mathrm{mod}\,q``.

We use the same symbol for both deterministic and random flattening, because the choice of either does not affect the algorithms that use them (other than making their results deterministic or random, respectively).


#### Polynomial flattening

*Paper:* Section 2.5.

*Implementation:* [`flatten()`](@ref SGFHE.flatten_poly)

When the flattening operation is applied to a polynomial ``a(x) \in R_{n,q}``, we decompose every coefficient separately and use the results to construct ``\ell`` polynomials:
```math
a(x) \triangleleft (B, \ell) = (u_0(x), \dots, u_{\ell-1}(x)) \in R_{n,q}^\ell,
```
where
```math
\mathcal{C}[u_j(x)]_i = \left( \mathcal{C}[a(x)]_i \triangleleft (B, \ell) \right)_j,
\quad i = 0, \dots, n-1,\,j = 0, \dots, \ell-1.
```


### [External product (``\odot``)](@id External-product)

*Paper:* Section 2.5.

*Implementation:* [`external_product()`](@ref SGFHE.external_product)

External product acts on a two-element vector of polynomials ``(a(x), b(x)) \in R_{n,q}^2`` and a matrix of polynomials ``A \in R_{n,q}^{2\ell \times 2}``:
```math
\boldsymbol{v} \odot (A, B, \ell)
  = \left(
    a(x) \triangleleft (B, \ell), b(x) \triangleleft (B, \ell)
  \right) A \in R_{n,q}^2,
```
where the results of applying ``\triangleleft`` to ``a(x)`` and ``b(x)`` (each an ``\ell``-vector of polynomials) are concatenated, producing a ``2\ell``-vector.


### [Shortened external product (``\odot``)](@id Shortened-external-product)

*Paper:* Section 2.5.

*Implementation:* [`shortened_external_product()`](@ref SGFHE.shortened_external_product)

This is a version of the external product acting on a single polynomial ``a(x) \in R_{n,q}`` and a matrix of polynomials ``A \in R_{n,q}^{2\ell \times 2}``:
```math
\boldsymbol{a} \odot (A, B, \ell) = \left( a(x) \triangleleft (B, \ell) \right) A_{\ell+1 \dots 2\ell, \dots} \in R_{n,q}^2,
```
where ``A_{\ell+1 \dots 2\ell, \dots}`` denotes the bottom half of the matrix ``A`` (the rows from ``\ell+1`` to ``2\ell``).


## [Scheme parameters](@id Scheme-parameters-theory)

*Paper:* Section 4.1.

*Implementation:* [`Params`](@ref Params)

The constants introduced here will be used throughout the rest of the section.

The main parameter of the FHE scheme is the polynomial length ``n \ge 64``, which must be a power of ``2``.

Other parameters that need to be chosen:

* LWE modulus is set to ``r = 16n`` (it is possible to pick any ``r >= 16n`` as long as ``r`` is divisible by ``8``, which will require several minor adjustments in the algorithms; see the paper for details).
* RLWE modulus ``q >= nr``. In this implementation we choose a prime ``q`` such that ``q-1`` is divisible by ``2n`` (to be able to multiply polynomials of length ``n`` using NTT).
* Decomposition length ``\ell = 2``. See [Flattening](@ref) for details. Note that this value is "hardcoded" in the bootstrap modulus ``Q`` below (``Q \le B^\ell``), so it cannot be trivially adjusted.
* Bootstrap modulus ``1220 r^4 n^2 \le Q \le 1225 r^4 n^2 = B^2`` (see below for the value of ``B``). Again, we choose a prime ``Q`` such that ``Q-1`` is divisible by ``2m`` (to be able to multiply polynomials of length ``m`` using NTT).

Dependent parameters:

* Extended polynomial length for bootstrapping ``m = r / 2``.
* Decomposition base ``B = 35 r^2 n``. See [Flattening](@ref) for details.
* LWE error level ``D_r = r/4``.
* RLWE error level ``D_q = q \div 4``.
* Bootstrap error level ``\tilde{D}_Q = Q \div 8``.
* Lower limit for the bit size of LWE elements ``t = \lceil \log_2 r \rceil - 1`` (so ``2^t < r \le 2^{t+1}``).

We will also need the gadget matrix ``G \in \mathbb{Z}_Q^{2\ell \times 2}``
```math
G = \begin{bmatrix}
    1 & 0 \\
    B & 0 \\
    \dots & 0 \\
    B^{\ell-1} & 0 \\
    0 & 1 \\
    0 & B \\
    0 & \dots \\
    0 & B^{\ell-1} \\
\end{bmatrix}
```


## Private key generation

*Paper:* Section 4.1.

*Implementation:* [`PrivateKey`](@ref PrivateKey)

A private key is simply a vector ``\boldsymbol{s} \in \{ 0, 1 \}^n`` with the elements ``s_i = \mathcal{U}(0, 1)``, ``i \in [0, n-1]``.


## Private key encryption to an RLWE

*Paper:* Fig. 2.

*Implementation:* [`encrypt`](@ref encrypt)

Given:

* A private key ``\boldsymbol{s} \in \{ 0, 1 \}^n``;
* An ``n``-bit message to encrypt ``\boldsymbol{m}  \in \{ 0, 1 \}^n``.

Algorithm:

* Generate ``\boldsymbol{u} \in \{ 0, 1 \}^n``, ``u_i = \mathcal{U}(0, 1)``.
* Expand it into a random polynomial:
  ```math
  \boldsymbol{a} = \mathrm{Expand}[\boldsymbol{u}, t+1], \\
  a(x) = \mathcal{P}[\boldsymbol{a}] \in R_{n,r}.
  ```
* Generate
  ```math
  \boldsymbol{w} \in \mathbb{Z}_r^n: \ w_i = \mathcal{U}(-D_r \div 8, D_r \div 8) \,\mathrm{mod}\,r, \\
  w(x) = \mathcal{P}[\boldsymbol{w}] \in R_{n,r}.
  ```
* Calculate
  ```math
  b(x) = a(x) s(x) + w(x) + \mathcal{P}[\boldsymbol{m}] D_r \,\mathrm{mod}\,(x^n+1,r).
  ```

Result: an RLWE cipher ``(a(x),b(x)) \in R_{n,r}^2``.


### [Space-optimal representation](@id Space-optimal-representation-private)

*Paper:* Fig. 2.

*Implementation:* [`encrypt_optimal`](@ref encrypt_optimal)

Starting from where the algorithm in the previous section ended,

* Calculate ``\tilde{b}(x) = b(x) \div 2^{t-4}`` (in other words, only highest 5 bits of ``b(x)``'s coefficients are important)
* Convert the vector ``\tilde{\boldsymbol{b}} = \mathcal{C}[\tilde{b}(x)]`` into a matrix ``V \in \{0, 1\}^{n \times 5}`` where the ``i``-th row is the bit representation of ``\tilde{b}_i``.

Result: a pair ``(\boldsymbol{u}, V) \in \{0, 1\}^n \times \{0, 1\}^{n \times 5}``.


### Restoring from the space-optimal representation

*Paper:* Fig. 2.

*Implementation:* [`normalize_ciphertext`](@ref normalize_ciphertext)

Given a pair ``(\boldsymbol{u}, V) \in \{0, 1\}^n \times \{0, 1\}^{n \times 5}``:

* Expand ``\boldsymbol{u}`` into a polynomial:
  ```math
  \boldsymbol{a} = \mathrm{Expand}[\boldsymbol{u}, t+1], \\
  a(x) = \mathcal{P}[\boldsymbol{a}] \in R_{n,r}.
  ```
* Convert ``V`` into a polynomial:
  ```math
  \tilde{\boldsymbol{b}} \in \mathbb{Z}_r^n: \ \tilde{b}_i = \sum_{j=0}^4 V_{ij} 2^j, \\
  b(x) = 2^{t-4} \mathcal{P}[\tilde{\boldsymbol{b}}]
  ```

Result: an RLWE cipher ``(a(x),b(x)) \in R_{n,r}^2``.


## Public key generation

*Paper:* Section 4.1.

*Implementation:* [`PrivateKey`](@ref PrivateKey)

Given a private key ``\boldsymbol{s} \in \{ 0, 1 \}^n``:

* Generate
  ```math
  \boldsymbol{k}^{(0)} \in \mathbb{Z}_q^n:\ k^{(0)}_i = \mathcal{U}(0, q-1).
  ```
* Generate
  ```math
  \boldsymbol{e} \in \mathbb{Z}_q^n:\ e_i = \mathcal{U}(-c, c),
  ```
  where ``c`` is the largest integer such that ``c < D_q / (41n)``.
* Calculate
  ```math
  k^{(1)}(x) = k^{(0)}(x) s(x) + e(x) \,\mathrm{mod}\,(x^n+1, q).
  ```

Result: a pair ``(k^{(0)}(x), k^{(1)}(x)) \in R_{n,q}^2``.


## Public key encryption to an RLWE

*Paper:* Fig. 3.

*Implementation:* [`encrypt()`](@ref encrypt)

Given:

* A public key ``(k^{(0)}(x), k^{(1)}(x)) \in R_{n,q}^2``;
* A message to encrypt ``\boldsymbol{m} \in \{0, 1\}^n``.

Algorithm:

* Generate
  ```math
  \boldsymbol{u} \in \mathbb{Z}_q^n:\ \mathcal{C}[u(x)]_i = \mathcal{U}(-1, 1), \\
  u(x) = \mathcal{P}[\boldsymbol{u}] \in R_{n,q}.
  ```
* Generate
  ```math
  \boldsymbol{w}^{(1)} \in \mathbb{Z}_q^n:\ w^{(1)}_i = \mathcal{U}(-D_q/(41n), D_q/(41n)) \,\mathrm{mod}\,q, \\
  w^{(1)}(x) = \mathcal{P}[\boldsymbol{w}^{(1)}] \in R_{n,q}.
  ```
* Generate
  ```math
  \boldsymbol{w}^{(2)} \in \mathbb{Z}_q^n:\ \mathcal{U}(-D_q/82, D_q/82) \,\mathrm{mod}\,q, \\
  w^{(2)}(x) = \mathcal{P}[\boldsymbol{w}^{(2)}] \in R_{n,q}.
  ```
* Calculate
  ```math
  a_1(x) = k^{(0)}(x) u(x) + w^{(1)}(x) \,\mathrm{mod}\,(x^n+1, q), \\
  b_1(x) = k^{(1)}(x) u(x) + w^{(2)}(x) + m(x) D_q \,\mathrm{mod}\,(x^n+1, q).
  ```
* Calculate
  ```math
  a(x) = \mathrm{ModRed}[a_1(x), q, r], \\
  b(x) = \mathrm{ModRed}[b_1(x), q, r].
  ```

Result: an RLWE cipher ``(a(x),b(x)) \in R_{n,r}^2``.


### Space-optimal representation

*Paper:* Fig. 3.

*Implementation:* [`encrypt_optimal()`](@ref encrypt_optimal)

Starting from where the algorithm in the previous section ended,

* Calculate ``\tilde{b}(x) = b(x) \div 2^{t-5}`` (in other words, only highest 6 bits of ``b(x)``'s coefficients are important)
* Convert the vector ``\tilde{\boldsymbol{b}} = \mathcal{C}[\tilde{b}(x)]`` into a matrix ``V \in \{0, 1\}^{n \times 6}`` where the ``i``-th row is the bit representation of ``\tilde{b}_i``.
* Convert the vector ``\boldsymbol{a} = \mathcal{C}[a(x)]`` into a matrix ``U \in \{0, 1\}^{n \times (t+1)}`` where the ``i``-th row is the bit representation of ``a_i``.

Result: a pair ``(U, V) \in \{0, 1\}^{n \times (t+1)} \times \{0, 1\}^{n \times 6}``.


### Restoring from the space-optimal representation

*Paper:* Fig. 3.

*Implementation:* [`normalize_ciphertext()`](@ref normalize_ciphertext)

Given a pair ``(U, V) \in \{0, 1\}^{n \times (t+1)} \times \{0, 1\}^{n \times 6}``:

* Convert ``U`` into a polynomial:
  ```math
  \boldsymbol{a} \in \mathbb{Z}_r^n:\ a_i = \sum_{j=0}^t U_{ij} 2^j, \\
  a(x) = \mathcal{P}[\boldsymbol{a}] \in R_{n,r}.
  ```
* Convert ``V`` into a polynomial:
  ```math
  \tilde{\boldsymbol{b}} \in \mathbb{Z}_r^n:\ \tilde{b}_i = \sum_{j=0}^t V_{ij} 2^j, \\
  b(x) = 2^{t-5} \mathcal{P}[\tilde{\boldsymbol{a}}] \in R_{n,r}.
  ```

Result: an RLWE cipher ``(a(x),b(x)) \in R_{n,r}^2``.


## Creating LWEs


### Private key encryption to an LWE

*Note:* this is not explicitly described in the paper, but Lemma 2.3 sets the maximum error one can introduce during encryption.

*Implementation:* [`encrypt()`](@ref encrypt)

Given:

* a message to encrypt ``y \in \{0, 1\}``;
* a secret key ``\boldsymbol{s} \in \{0, 1\}^n``.

Algorithm:

* Generate a vector ``\boldsymbol{a}^\prime \in \mathbb{Z}_q^n`` with the elements ``a^\prime_i = \mathcal{U}(0, q-1)``.
* Generate ``e = \mathcal{U}(-\tau, \tau) \,\mathrm{mod}\, q`` where ``\tau = (q (n-3)) \div (2r)``.
* Calculate
  ```math
  b^\prime = \boldsymbol{a} \cdot \boldsymbol{s} + e + y D_q \,\mathrm{mod}\, q.
  ```
* Modulus reduction:
  ```math
  \boldsymbol{a} = \mathrm{ModRed}[\boldsymbol{a}^\prime, q, r], \\
  b = \mathrm{ModRed}[b^\prime, q, r].
  ```

Result: an LWE cipher ``\mathrm{LWE}_{\boldsymbol{s}}(y) = (\boldsymbol{a}, b) \in \mathbb{Z}_r^n \times \mathbb{Z}_r``.


### Extracting an LWE from a packed RLWE

*Paper:* Section 4.5.1.

*Implementation:* [`split_ciphertext()`](@ref split_ciphertext)

Given an RLWE cipher ``\mathrm{RLWE}_{\boldsymbol{s}}(\boldsymbol{y}) = (a(x), b(x)) \in R_{n,r}^2`` encrypting an ``n``-bit message ``\boldsymbol{y} \in \{0, 1\}^n``, the LWE ``\mathrm{LWE}_{\boldsymbol{s}}(y_i)`` encrpyting the ``i``-th bit is
```math
\mathrm{LWE}_{\boldsymbol{s}}(y_i) = (\mathrm{Extract}[\mathcal{C}[a(x)], i, n], \mathcal{C}[b(x)]_i).
```


### Extracting an LWE from a RLWE

*Note:* this is not described in the original paper.

*Implementation:* [`split_ciphertext()`](@ref split_ciphertext)

The procedure is the same as for the packed RLWE. Given an RLWE cipher ``\mathrm{RLWE}_{\boldsymbol{s}}(\boldsymbol{y}) = (a(x), b(x)) \in R_{m,r}^2`` encrypting an ``n``-bit message ``\boldsymbol{y} \in \{0, 1\}^n``, the LWE ``\mathrm{LWE}_{\boldsymbol{s}}(y_i)`` encrpyting the ``i``-th bit is
```math
\mathrm{LWE}_{\boldsymbol{s}}(y_i) = (\mathrm{Extract}[\mathcal{C}[a(x)], i, n], \mathcal{C}[b(x)]_i).
```


## Decryption


### Decrypting an LWE

*Paper:* Section 2.3.

*Implementation:* [`decrypt()`](@ref decrypt)

Given:

* An LWE cipher encrypting ``y \in \{0,1\}`` ``\mathrm{LWE}_{\boldsymbol{s}}(y) = (\boldsymbol{a}, b) \in \mathbb{Z}_r^n \times \mathbb{Z}_r``;
* A secret key ``\boldsymbol{s} \in \{0, 1\}^n``.

Algorithm:

* Calculate
  ```math
  b_1 = b - \boldsymbol{a} \cdot \boldsymbol{s} + D_r / 2 \,\mathrm{mod}\,r.
  ```
* Calculate ``y = b_1 \div D_r``.

The addition of ``D_r / 2`` as compared to the algorithm provided in the paper allows us to use ``\div`` (that is, the regular integer floor division by ``D_r``) as compared to ``\lfloor \rceil`` in the paper.


### Decrypting a packed RLWE

*Paper:* Section 4.2.

*Implementation:* [`decrypt()`](@ref decrypt)

Given:

* An RLWE cipher encrypting ``\boldsymbol{y} \in \{0,1\}^n`` ``\mathrm{RLWE}_{\boldsymbol{s}}(\boldsymbol{y}) = (a(x), b(x)) \in R_{n,r}^2``;
* A secret key ``\boldsymbol{s} \in \{0, 1\}^n``.

Algorithm:

* Calculate
  ```math
  b_1(x) = b(x) - a(x) s(x) \,\mathrm{mod}\,(x^n+1, r).
  ```
* Calculate
  ```math
  y_i = (\mathcal{C}[b_1(x)]_i + D_r / 2 \,\mathrm{mod}\,r) \div D_r,
  \quad i = 0, \dots, n-1.
  ```

Again, the shift by ``D_r / 2`` allows us to use the floor division instead of rounding division.


### Decrypting an RLWE

*Note:* this algorithm is not given in the original paper.

*Implementation:* [`decrypt()`](@ref decrypt)

The algorithm is almost identical to the one for a packed RLWE, with the only change being the polynomial modulus in the first step.

Given:

* An RLWE cipher encrypting ``\boldsymbol{y} \in \{0,1\}^n`` ``\mathrm{RLWE}_{\boldsymbol{s}}(\boldsymbol{y}) = (a(x), b(x)) \in R_{m,r}^2``;
* A secret key ``\boldsymbol{s} \in \{0, 1\}^n``.

Algorithm:

* Calculate
  ```math
  b_1(x) = b(x) - a(x) s(x) \,\mathrm{mod}\,(x^m+1, r).
  ```
* Calculate
  ```math
  y_i = (\mathcal{C}[b_1(x)]_i + D_r / 2 \,\mathrm{mod}\,r) \div D_r,
  \quad i = 0, \dots, n-1.
  ```


## Bootstrap key generation

*Paper:* Section 4.1.

*Implementation:* [`BootstrapKey`](@ref BootstrapKey)

Given a secret key ``\boldsymbol{s} \in \{0, 1\}^n``, for each ``i = 0, \dots, n-1``:

* Generate a matrix of polynomials ``A \in R_{m,Q}^{2\ell \times n}`` where
  ```math
  \mathcal{C}[A_{ji}(x)]_k = \mathcal{U}(0, Q-1), \quad k = 0, \dots, n-1.
  ```
* Generate a matrix of polynomials ``E \in R_{m,Q}^{2\ell \times n}`` where
  ```math
  \mathcal{C}[E_{ji}(x)]_k = \mathcal{U}(-n, n) \,\mathrm{mod}\,Q, \quad k = 0, \dots, n-1.
  ```
* Calculate
  ```math
  B(x) = A(x) \mathcal{P}[\boldsymbol{s}](x) + E(x) \,\mathrm{mod}\,(x^m+1, Q).
  ```
* For ``i = 0, \dots, n-1``, calculate
  ```math
  C^{(i)} = \begin{bmatrix}
      A_{1i}(x) & B_{1i}(x) \\
      A_{2i}(x) & B_{2i}(x) \\
      A_{3i}(x) & B_{3i}(x) \\
      A_{4i}(x) & B_{4i}(x)
  \end{bmatrix} + s_i G \,\mathrm{mod}\,Q.
  ```
  where ``G`` is the gadget matrix (see [Scheme parameters](@ref Scheme-parameters-theory)).

The resulting bootstrap key is a list of ``n`` matrices of polynomials ``C^{(i)} \in R_{m,Q}^{2\ell \times n}``.


## [Bootstrap](@id Bootstrap-theory)

*Paper:* Fig. 1.

*Implementation:* [`bootstrap()`](@ref bootstrap)

Given:

* A bootstrap key --- a list of ``n`` matrices of polynomials ``C^{(i)} \in R_{m,Q}^{2\ell \times n}``;
* Two LWE ciphers ``\boldsymbol{v}_1 = \mathrm{LWE}_{\boldsymbol{s}}(y_1)``, ``\boldsymbol{v}_2 = \mathrm{LWE}_{\boldsymbol{s}}(y_2)``, where ``\boldsymbol{v}_1, \boldsymbol{v}_2 \in \mathbb{Z}^n_r \times \mathbb{Z}_r``.

Algorithm:

* Sum both LWE ciphers element-wise:
  ```math
  \boldsymbol{u} = \boldsymbol{v}_1 + \boldsymbol{v}_2 = (u_0, \dots, u_{n-1}, u_n) \in \mathbb{Z}^n_r \times \mathbb{Z}_r.
  ```
* Build
  ```math
  t(x) = \sum_{j=-D_r+1}^{D_r-1} x^j \,\mathrm{mod}\,(x^m+1, Q).
  ```
  See [Negative powers in polynomials](@ref) for the explanation of how the modulus is taken.
* Initilize ``a(x), b(x) \in R_{m,Q}``
  ```math
  a(x) = 0, \\
  b(x) = t(x) x^{-u_n} \tilde{D}_Q \,\mathrm{mod}\,(x^m+1, Q).
  ```
* For ``i = 0 \dots n-1``:
  ```math
  A^{(i)} = G + (x^{u_i} - 1) C^{(i)} \,\mathrm{mod}\,(x^m+1, Q), \\
  (a(x), b(x)) \leftarrow (a(x), b(x)) \odot (A^{(i)}, B, \ell),
  ```
  where G is the gadget matrix (see [Scheme parameters](@ref Scheme-parameters-theory)) and ``\odot`` is the [external product](@ref External-product).
* Build LWEs
  ```math
  \boldsymbol{a}^{\mathrm{AND}} = (\mathrm{Extract}[\boldsymbol{a}, 3m/4], \tilde{D}_Q + b_{3m/4}), \\
  \boldsymbol{a}^{\mathrm{OR}} = (-\mathrm{Extract}[\boldsymbol{a}, 3m/4], \tilde{D}_Q - b_{m/4}), \\
  \boldsymbol{a}^{\mathrm{XOR}} = \boldsymbol{a}_{\mathrm{OR}} - \boldsymbol{a}_{\mathrm{AND}}.
  ```
  where
  ```math
  \boldsymbol{a} = \mathcal{C}[a(x)] \in \mathbb{Z}_Q^m, \\
  \boldsymbol{b} = \mathcal{C}[b(x)] \in \mathbb{Z}_Q^m.
  ```
* Modulus reduction:
  ```math
  \boldsymbol{c}^{\mathrm{AND}} = \mathrm{ModRed}[\boldsymbol{a}^{\mathrm{AND}}, Q, r], \\
  \boldsymbol{c}^{\mathrm{OR}} = \mathrm{ModRed}[\boldsymbol{a}^{\mathrm{OR}}, Q, r], \\
  \boldsymbol{c}^{\mathrm{XOR}} = \mathrm{ModRed}[\boldsymbol{a}^{\mathrm{XOR}}, Q, r].
  ```
  Here ``\boldsymbol{c}^{\mathrm{AND}}, \boldsymbol{c}^{\mathrm{OR}}, \boldsymbol{c}^{\mathrm{XOR}} \in \mathbb{Z}_r^n \times \mathbb{Z}_r``.

Result: LWE ciphers of ``y_1 \wedge y_2`` (AND), ``y_1 \vee y_2`` (OR) and ``y_1 \oplus y_2`` (XOR):
```math
\boldsymbol{c}^{\mathrm{AND}} = \mathrm{LWE}_{\boldsymbol{s}}(y_1 \wedge y_2), \\
\boldsymbol{c}^{\mathrm{OR}} = \mathrm{LWE}_{\boldsymbol{s}}(y_1 \vee y_2), \\
\boldsymbol{c}^{\mathrm{XOR}} = \mathrm{LWE}_{\boldsymbol{s}}(y_1 \oplus y_2).
```

One can use either the deterministic or the random version of the external product ``\odot``, depending on the desired behavior of the bootstrap.


## Packing LWEs

*Paper:* Lemma 4.3.

*Implementation:* [`pack_encrypted_bits()`](@ref pack_encrypted_bits)

Given:
* ``n`` LWE ciphers ``\boldsymbol{z}^{(i)} = \mathrm{LWE}_{\boldsymbol{s}}(y_i) \in \mathbb{Z}_r^n \times \mathbb{Z}_r, i = 0, \dots, n-1``;
* a bootstrap key ``\{C^{(i)}\}, i = 0, \dots, n-1``.

Algorithm:

* Create a trivial LWE cipher encrypting ``1``: ``\mathrm{LWE}_{\boldsymbol{s}}(1) = (\boldsymbol{0}, D_r) \in \mathbb{Z}_r^n \times \mathbb{Z}_r``.
* For ``i = 0, \dots, n-1`` apply the [Bootstrap](@ref) algorithm to the pair ``E_{\boldsymbol{s}}(1)`` and ``\boldsymbol{z}^{(i)}`` up to and not including the modulus reduction step, leaving only the first result (corresponding to the AND gate). This will produce ``n`` LWEs:
  ```math
  \boldsymbol{a}^{\mathrm{AND},(i)}
  \equiv (\boldsymbol{a}^{(i)}, b^{(i)})
  \in \mathbb{Z}_Q^n \times \mathbb{Z}_Q.
  ```
* For ``i = 0, \dots, n-1``, build polynomials
  ```math
  \tilde{a}^{(i)}(x) = \mathcal{P}[\tilde{\boldsymbol{a}}^{(i)}] \in R_{m,Q},
  ```
  where ``\tilde{a}^{(i)}_j = a^{(j)}_i``, ``j = 0, \dots, n-1``. Note that we are treating the resulting polynomials as modulo ``(x^m+1)`` despite them only having powers of ``x`` up to ``n-1``.
* Create a polynomial ``\tilde{b}(x) = \mathcal{P}[\boldsymbol{b}] \in R_{m,Q}``. Again, we are treating the resulting polynomial as modulo ``(x^m+1)``.
* Calculate:
  ```math
  \left( \tilde{w}(x), \tilde{v}(x) \right) = \sum_{i=0}^{n-1}
    a^{(i)}(x) \odot (C^{(i)}, B, \ell),
  ```
  where ``\odot`` is the [shortened external product](@ref Shortened-external-product).
* Calculate:
  ```math
  \tilde{w}_1(x) = -\tilde{w}(x), \\
  \tilde{v}_1(x) = b(x) - \tilde{v}(x).
  ```
* Modulus reduction:
  ```math
  w(x) = \mathrm{ModRed}[\tilde{w}_1(x), Q, r], \\
  v(x) = \mathrm{ModRed}[\tilde{v}_1(x), Q, r].
  ```

Result: an RLWE cipher ``(w(x), v(x)) \in R_{m,r}^2`` encrypting the vector ``\boldsymbol{y}``.
