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
        return defaultMDConfig;
    }

    // eslint-disable-next-line class-methods-use-this
    get uiSchema() {
        return {
            nstep: {},
            dt: {},
            emass: {},
            tempw: {},
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
                    title: "Effective electron mass",
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
