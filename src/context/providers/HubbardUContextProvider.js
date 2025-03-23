import {
    JSONSchemaFormDataProvider,
    MaterialContextMixin,
    MethodDataContextMixin,
} from "@mat3ra/code/dist/js/context";
import { Made } from "@mat3ra/made";
import { mix } from "mixwith";

import { sortArrayByOrder } from "../../utils";

const defaultHubbardConfig = {
    atomicSpecies: "",
    atomicOrbital: "3d",
    hubbardUValue: 1.0,
};

export class HubbardUContextProvider extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
    MethodDataContextMixin,
) {
    static Material = Made.Material;

    constructor(config) {
        super(config);
        this.uniqueElements = this.material?.Basis?.uniqueElements || [];
        // orbitals are sorted according to stability (Madelung's rule)
        this.orbitalList = [
            "1s",
            "2s",
            "2p",
            "3s",
            "3p",
            "4s",
            "3d",
            "4p",
            "5s",
            "4d",
            "5p",
            "6s",
            "4f",
            "5d",
            "6p",
            "7s",
            "5f",
            "6d",
            "7p",
            "8s",
        ];
        const _elementsWithLabels = this.material?.Basis?.elementsWithLabelsArray || [];
        this.uniqueElementsWithLabels = [...new Set(_elementsWithLabels)];
        this.firstElement =
            this.uniqueElementsWithLabels?.length > 0 ? this.uniqueElementsWithLabels[0] : "";
    }

    get defaultData() {
        const valenceOrbitals = this.getValenceOrbitals(this.firstElement);
        return [
            {
                ...defaultHubbardConfig,
                atomicSpecies: this.firstElement,
                atomicOrbital:
                    valenceOrbitals.length > 0
                        ? valenceOrbitals[valenceOrbitals.length - 1]
                        : defaultHubbardConfig.atomicOrbital,
            },
        ];
    }

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: false,
                removable: true,
            },
            items: {
                atomicSpecies: this.defaultFieldStyles,
                atomicOrbital: this.defaultFieldStyles,
                hubbardUValue: this.defaultFieldStyles,
            },
        };
    }

    getValenceOrbitals = (element) => {
        const pseudos = this.methodData?.pseudo || [];
        let valenceConfig = [];
        pseudos.every((data) => {
            if (data.element === element) {
                valenceConfig = data?.valenceConfiguration || [];
            }
            return data.element !== element; // break when first match is found
        });
        const valenceOrbitals = valenceConfig.map((item) => item.orbitalName.toLowerCase());
        return sortArrayByOrder(valenceOrbitals, this.orbitalList);
    };

    getElementSymbol = (elementWithLabel) => {
        // exclude single digit label in the end of symbol if present
        // 1 is added to the label below to take care of possible label 0
        return parseInt(elementWithLabel.slice(-1), 10) + 1
            ? elementWithLabel.slice(0, -1)
            : elementWithLabel;
    };

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            title: "",
            description: "Hubbard U parameters for DFT+U or DFT+U+V calculation.",
            type: "array",
            items: {
                type: "object",
                properties: {
                    atomicSpecies: {
                        type: "string",
                        title: "Atomic species",
                        enum: this.uniqueElementsWithLabels,
                        default: this.firstElement,
                    },
                    atomicOrbital: {
                        type: "string",
                        title: "Atomic orbital",
                    },
                    hubbardUValue: {
                        type: "number",
                        title: "Hubbard U (eV)",
                        default: defaultHubbardConfig.hubbardUValue,
                    },
                },
                dependencies: {
                    atomicSpecies: {
                        oneOf: this.uniqueElementsWithLabels.map((elementWithLabel) => {
                            const orbitals = this.getValenceOrbitals(
                                this.getElementSymbol(elementWithLabel),
                            );
                            return {
                                properties: {
                                    atomicSpecies: {
                                        enum: [elementWithLabel],
                                    },
                                    atomicOrbital: {
                                        enum: orbitals.length > 0 ? orbitals : this.orbitalList,
                                        default:
                                            orbitals.length > 0
                                                ? orbitals[orbitals.length - 1]
                                                : defaultHubbardConfig.atomicOrbital,
                                    },
                                },
                            };
                        }),
                    },
                },
            },
        };
    }
}
