# dangerous

Small Julia library to simulate 1D and 2D nuclear magnetic resonance pulse programmes.


> [!WARNING]
> This is a toy project and is specifically named to reflect that.
> Please don't take it seriously.
> In particular, please don't use it.

## Why is it dangerous?

I ran the phrase 'nuclear magnetic resonance' through Google Translate 500 times and got back 'it is dangerous'.


## What can you do with it?

If you disobey the instructions in the warning above, you can sometimes generate nice Lorentzian lineshapes.

```shell
cd examples
julia --project .
]instantiate
```

then

```julia
julia> include("single_spin.jl")
```

![1 spin](https://github.com/user-attachments/assets/473a3535-c4fe-47d8-a2b4-6d0bd25588b2)


Alternatively, if you want a pair of weakly-coupled spins (i.e. an AX spin system):

```julia
julia> include("ax.jl")
```

![2 spins](https://github.com/user-attachments/assets/82a35c7a-fcf2-4c23-ad7b-d9148ad5c1f4)

