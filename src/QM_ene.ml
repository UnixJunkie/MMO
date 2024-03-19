module QM_ene : sig
  type t

  val of_pyobject : Pytypes.pyobject -> t
  val to_pyobject : t -> Pytypes.pyobject
  val __init__ : anums:int array -> unit -> t

  val score_conf :
    t -> xs:float array -> ys:float array -> zs:float array -> unit -> float
end = struct
  let filter_opt l = List.filter_map Fun.id l

  let py_module =
    lazy
      (let source =
         {pyml_bindgen_string_literal|# Python interface; will be wrapped by an OCaml module
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
|pyml_bindgen_string_literal}
       in
       let filename =
         {pyml_bindgen_string_literal|QM_ene.py|pyml_bindgen_string_literal}
       in
       let bytecode = Py.compile ~filename ~source `Exec in
       Py.Import.exec_code_module
         {pyml_bindgen_string_literal|QM_ene|pyml_bindgen_string_literal}
         bytecode)

  let import_module () = Lazy.force py_module

  type t = Pytypes.pyobject

  let of_pyobject pyo = pyo
  let to_pyobject x = x

  let __init__ ~anums () =
    let callable = Py.Module.get (import_module ()) "Torchani" in
    let kwargs =
      filter_opt [ Some ("anums", Py.List.of_array_map Py.Int.of_int anums) ]
    in
    of_pyobject @@ Py.Callable.to_function_with_keywords callable [||] kwargs

  let score_conf t ~xs ~ys ~zs () =
    let callable = Py.Object.find_attr_string t "score_conf" in
    let kwargs =
      filter_opt
        [
          Some ("xs", Py.List.of_array_map Py.Float.of_float xs);
          Some ("ys", Py.List.of_array_map Py.Float.of_float ys);
          Some ("zs", Py.List.of_array_map Py.Float.of_float zs);
        ]
    in
    Py.Float.to_float
    @@ Py.Callable.to_function_with_keywords callable [||] kwargs
end
