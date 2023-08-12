import { Application } from "@exabyte-io/ade.js";
import { ApplicationContextMixin, ContextProvider } from "@exabyte-io/code.js/dist/context";
import { mix } from "mixwith";

const defaultHubbardConfig = {
    atomicSpecies: "",
    atomicOrbital: "",
    hubbardUValue: 0.01,
};

export class HubbardContextProvider extends mix(ContextProvider).with(ApplicationContextMixin) {
    static Application = Application;

    // eslint-disable-next-line class-methods-use-this
    get uiSchema() {
        return {
            atomicSpecies: {},
            atomicOrbital: {},
            hubbardUValue: {},
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return defaultHubbardConfig;
    }

    // eslint-disable-next-line class-methods-use-this
    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U (Quantum Espresso) calculation.",
            type: "object",
            properties: {
                atomicSpecies: {
                    type: "string",
                    default: defaultHubbardConfig.atomicSpecies,
                },
                atomicOrbital: {
                    type: "string",
                    default: defaultHubbardConfig.atomicOrbital,
                },
                hubbardUValue: {
                    type: "number",
                    default: defaultHubbardConfig.hubbardUValue,
                },
            },
        };
    }
}
