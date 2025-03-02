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
    atomicOrbital: "2p",
    hubbardUValue: 1.0,
};

// Madelung's rule
const orbitalsByStability = [
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

export class HubbardUContextProvider extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
    MethodDataContextMixin,
) {
    static Material = Made.Material;

    constructor(config) {
        super(config);
        this.uniqueElements = this.material?.Basis?.uniqueElements || [];
        this.orbitalList = [
            "2p",
            "3s",
            "3p",
            "3d",
            "4s",
            "4p",
            "4d",
            "4f",
            "5s",
            "5p",
            "5d",
            "5f",
            "6s",
            "6p",
            "6d",
            "7s",
            "7p",
            "7d",
        ];
        const _elementsWithLabels = this.material?.Basis?.elementsWithLabelsArray || [];
        this.uniqueElementsWithLabels = [...new Set(_elementsWithLabels)];
        this.firstElement =
            this.uniqueElementsWithLabels?.length > 0 ? this.uniqueElementsWithLabels[0] : "";
    }

    get defaultData() {
        return [
            {
                ...defaultHubbardConfig,
                atomicSpecies: this.firstElement,
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

    _getValenceOrbitals = (element) => {
        const pseudos = this.methodData?.pseudo || [];
        let valenceConfig;
        pseudos.forEach((data) => {
            if (data.element === element) {
                valenceConfig = data.valenceConfiguration || [];
            }
        });
        const valenceOrbitals = valenceConfig.map((item) => item.orbitalName.toLowerCase());
        return sortArrayByOrder(valenceOrbitals, orbitalsByStability);
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
                    },
                    hubbardUValue: {
                        type: "number",
                    },
                    dependencies: {
                        atomicSpecies: {
                            oneOf: this.uniqueElementsWithLabels.map((elementWithLabel) => {
                                const element = parseInt(elementWithLabel.slice(-1), 10)
                                    ? elementWithLabel.slice(0, -1)
                                    : elementWithLabel;
                                return {
                                    properties: {
                                        atomicSpecies: {
                                            enum: [elementWithLabel],
                                        },
                                        atomicOrbital: {
                                            type: "string",
                                            title: "Atomic orbital",
                                            enum:
                                                this._getValenceOrbitals(element).length > 0
                                                    ? this._getValenceOrbitals(element)
                                                    : this.orbitalList,
                                            default:
                                                this._getValenceOrbitals(element).length > 0
                                                    ? this._getValenceOrbitals(element)[
                                                          this._getValenceOrbitals(element).length -
                                                              1
                                                      ]
                                                    : defaultHubbardConfig.atomicOrbital,
                                        },
                                        hubbardUValue: {
                                            type: "number",
                                            title: "Hubbard U (eV)",
                                            default: defaultHubbardConfig.hubbardUValue,
                                        },
                                    },
                                };
                            }),
                        },
                    },
                },
            },
        };
    }
}
