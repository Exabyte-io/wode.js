import { ConvergenceParameter } from "./parameter";

export class NonUniformKGridConvergence extends ConvergenceParameter {
    get increment() {
        return `[ ${this.name}[i] + K_INCREMENT[i] for i in range(len(${this.name}))]`;
    }

    get unitContext() {
        return {
            kgrid: {
                dimensions: [`{{${this.name}[0]}}`, `{{${this.name}[1]}}`, `{{${this.name}[2]}}`],
                shifts: [0, 0, 0],
            },
            isKgridEdited: true,
            isUsingJinjaVariables: true,
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get subworkflowContext() {
        return {
            K_INCREMENT: [1, 1, 1],
        };
    }
}
