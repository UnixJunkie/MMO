# Python interface; will be wrapped by an OCaml module
# compute ligand internal energy using a DNN approximation of a QM FF

import sys, torch, torchani

cpu_or_gpu = torch.device('cpu') # FBR: faster on my computer
# cpu_or_gpu = None
# if torch.cuda.is_available():
#     print("I: running on GPU", file=sys.stderr)
#     cpu_or_gpu = torch.device('cuda')
# else:
#     print("W: running on CPU", file=sys.stderr)
#     cpu_or_gpu = torch.device('cpu')
model = torchani.models.ANI2x(periodic_table_index=True).to(cpu_or_gpu)
hartree_to_kcal_per_mol = 627.5094740631

# # coordinates(Angstrom) and chemical species(anum); energies in Hartree
# coordinates = torch.tensor(xyz, requires_grad=False, device=cpu_or_gpu)
# species = torch.tensor(anums, device=cpu_or_gpu)
# energy = hartree_to_kcal_per_mol * model((species, coordinates)).energies.item()

class Torchani:
    # this is needed because the OCaml side needs to know how
    # to get an object of type t
    def __init__(self, anums):
        self.species = torch.tensor([anums], device=cpu_or_gpu)

    # energy for given conformer
    def score_conf(self, xs, ys, zs) -> float:
        coords = []
        n = len(xs)
        for i in range(n):
            coords.append([xs[i], ys[i], zs[i]])
        coordinates = torch.tensor([coords], requires_grad=False,
                                   device=cpu_or_gpu)
        return (hartree_to_kcal_per_mol * \
                model((self.species, coordinates)).energies.item())
