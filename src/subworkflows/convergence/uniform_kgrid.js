import { ConvergenceParameter } from "./parameter";

export class UniformKGridConvergence extends ConvergenceParameter {
    get increment() {
        return `${this.name} + ${this._increment}`;
    }

    get unitContext() {
        return {
            kgrid: {
                dimensions: [`{{${this.name}}}`, `{{${this.name}}}`, `{{${this.name}}}`],
                shifts: [0, 0, 0],
            },
            isKgridEdited: true,
            isUsingJinjaVariables: true,
        };
    }

    get finalValue() {
        return `${this.name} + 0`;
    }
}
