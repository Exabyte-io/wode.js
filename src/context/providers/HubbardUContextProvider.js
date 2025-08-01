import JSONSchemaFormDataProvider from "@exabyte-io/ade.js/dist/js/context/JSONSchemaFormDataProvider";
import { Made } from "@mat3ra/made";
import { Utils } from "@mat3ra/utils";

import { materialContextMixin } from "../mixins/MaterialContextMixin";
import { methodDataContextMixin } from "../mixins/MethodDataContextMixin";

const defaultHubbardConfig = {
    atomicSpecies: "",
    atomicOrbital: "3d",
    hubbardUValue: 1.0,
};

export class HubbardUContextProvider extends JSONSchemaFormDataProvider {
    constructor(config) {
        super(config);

        this.initMaterialContextMixin();
        this.initMethodDataContextMixin();

        this.uniqueElements = this.material?.Basis?.uniqueElements || [];
        // orbitals are sorted according to stability (Madelung's rule)
        this.orbitalListByStability = [
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

    getValenceOrbitalsByElement = (element) => {
        const valenceOrbitals = this.valenceOrbitals || [];
        let orbitals = [];
        valenceOrbitals.every((item) => {
            if (item.element === element) {
                orbitals = item?.valenceOrbitals || [];
            }
            return item.element !== element; // break when first match is found
        });

        return Utils.array.sortArrayByOrder(orbitals, this.orbitalListByStability);
    };

    orbitalDependencyArray = (elementList, atomicSpecies, atomicOrbital) => {
        return {
            oneOf: elementList.map((elementWithLabel) => {
                const orbitals = this.getValenceOrbitalsByElement(
                    Made.Basis.stripLabelToGetElementSymbol(elementWithLabel),
                );
                return {
                    properties: {
                        [atomicSpecies]: {
                            enum: [elementWithLabel],
                        },
                        [atomicOrbital]: {
                            enum: orbitals.length > 0 ? orbitals : this.orbitalListByStability,
                            default: this.getOutermostOrbital(orbitals),
                        },
                    },
                };
            }),
        };
    };

    getOutermostOrbital = (orbitals, defaultOrbital = defaultHubbardConfig.atomicOrbital) => {
        return orbitals.length > 0 ? orbitals[orbitals.length - 1] : defaultOrbital;
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
                    atomicSpecies: this.orbitalDependencyArray(
                        this.uniqueElementsWithLabels,
                        "atomicSpecies",
                        "atomicOrbital",
                    ),
                },
            },
        };
    }
}

materialContextMixin(HubbardUContextProvider.prototype);
methodDataContextMixin(HubbardUContextProvider.prototype);
