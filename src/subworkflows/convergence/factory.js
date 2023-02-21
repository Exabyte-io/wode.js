import { NonUniformKGridConvergence } from "./non_uniform_kgrid";
import { ConvergenceParameter } from "./parameter";
import { UniformKGridConvergence } from "./uniform_kgrid";

export function createConvergenceParameter({ name, initialValue }) {
    switch (name) {
        case "N_k":
            return new UniformKGridConvergence({ name, initialValue });
        case "N_k_xyz":
            return new NonUniformKGridConvergence({ name, initialValue });
        default:
            return new ConvergenceParameter({ name, initialValue });
    }
}
