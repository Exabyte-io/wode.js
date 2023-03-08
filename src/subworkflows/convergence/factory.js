import { NonUniformKGridConvergence } from "./non_uniform_kgrid";
import { ConvergenceParameter } from "./parameter";
import { UniformKGridConvergence } from "./uniform_kgrid";

export function createConvergenceParameter({ name, initialValue, increment }) {
    switch (name) {
        case "N_k":
            return new UniformKGridConvergence({ name, initialValue, increment });
        case "N_k_nonuniform":
            return new NonUniformKGridConvergence({ name, initialValue, increment });
        default:
            return new ConvergenceParameter({ name, initialValue, increment });
    }
}
