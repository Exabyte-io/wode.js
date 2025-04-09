import { HubbardUContextProvider } from "./HubbardUContextProvider";

const defaultHubbardConfig = {
    paramType: "U",
    atomicSpecies: "",
    atomicOrbital: "3d",
    value: 1.0,
};

export class HubbardJContextProvider extends HubbardUContextProvider {
    get defaultData() {
        return [
            {
                ...defaultHubbardConfig,
                atomicSpecies: this.firstElement,
                atomicOrbital: this.getOutermostOrbital(
                    this.getValenceOrbitalsByElement(this.firstElement),
                ),
            },
        ];
    }

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: true,
                removable: true,
            },
            items: {
                paramType: this.defaultFieldStyles,
                atomicSpecies: this.defaultFieldStyles,
                atomicOrbital: this.defaultFieldStyles,
                value: this.defaultFieldStyles,
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U+J calculation.",
            type: "array",
            items: {
                type: "object",
                properties: {
                    paramType: {
                        type: "string",
                        title: "Species",
                        enum: ["U", "J", "B", "E2", "E3"],
                        default: defaultHubbardConfig.paramType,
                    },
                    atomicSpecies: {
                        type: "string",
                        title: "Species",
                        enum: this.uniqueElementsWithLabels,
                        default: this.firstElement,
                    },
                    atomicOrbital: {
                        type: "string",
                        title: "Orbital",
                    },
                    value: {
                        type: "number",
                        title: "Value (eV)",
                        default: defaultHubbardConfig.value,
                    },
                },
                dependencies: {
                    atomicSpecies: this.orbitalDependencyArray(
                        this.uniqueElementsWithLabels,
                        "atomicSpecies",
                        "atomicOrbital",
                    ),
                },
            },
            minItems: 1,
        };
    }
}
