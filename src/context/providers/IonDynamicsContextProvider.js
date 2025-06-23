import JSONSchemaFormDataProvider from "@exabyte-io/ade.js/dist/js/context/JSONSchemaFormDataProvider";

const defaultMDConfig = {
    numberOfSteps: 100,
    timeStep: 5.0,
    electronMass: 100.0,
    temperature: 300.0,
};

export class IonDynamicsContextProvider extends JSONSchemaFormDataProvider {
    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return defaultMDConfig;
    }

    // eslint-disable-next-line class-methods-use-this
    get uiSchema() {
        return {
            numberOfSteps: {},
            timeStep: {},
            electronMass: {},
            temperature: {},
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            type: "object",
            description: "Important parameters for molecular dynamics calculation",
            properties: {
                numberOfSteps: {
                    type: "integer",
                    title: "numberOfSteps",
                    default: defaultMDConfig.numberOfSteps,
                },
                timeStep: {
                    type: "number",
                    title: "timeStep (Hartree a.u.)",
                    default: defaultMDConfig.timeStep,
                },
                electronMass: {
                    type: "number",
                    title: "Effective electron mass",
                    default: defaultMDConfig.electronMass,
                },
                temperature: {
                    type: "number",
                    title: "Ionic temperature (K)",
                    default: defaultMDConfig.temperature,
                },
            },
        };
    }
}
