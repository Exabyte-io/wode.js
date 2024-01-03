import { JSONSchemaFormDataProvider } from "@exabyte-io/code.js/dist/context";

const defaultMDConfig = {
    nstep: 100,
    dt: 5.0,
    emass: 100.0,
    tempw: 300.0,
};

export class MolDynContextProvider extends JSONSchemaFormDataProvider {
    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return [defaultMDConfig];
    }

    get uiSchemaStyled() {
        return {
            title: {
                "ui:classNames": "col-xs-12",
            },
            items: {
                nstep: this.defaultFieldStyles,
                dt: this.defaultFieldStyles,
                emass: this.defaultFieldStyles,
                tempw: this.defaultFieldStyles,
            },
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            type: "object",
            description: "Important parameters for molecular dynamics calculation",
            properties: {
                nstep: {
                    type: "integer",
                    title: "nstep",
                    default: defaultMDConfig.nstep,
                },
                dt: {
                    type: "number",
                    title: "dt (Hartree a.u.)",
                    default: defaultMDConfig.dt,
                },
                emass: {
                    type: "number",
                    title: "emass (a.u.)",
                    description: "Effective electron mass in atomic units",
                    default: defaultMDConfig.emass,
                },
                tempw: {
                    type: "number",
                    title: "Ionic temperature (K)",
                    default: defaultMDConfig.tempw,
                },
            },
        };
    }
}
