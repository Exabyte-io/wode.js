import { HubbardUContextProvider } from "./HubbardUContextProvider";

const defaultHubbardConfig = {
    atomicSpecies: "",
    atomicOrbital: "3d",
    atomicSpecies2: "",
    atomicOrbital2: "3d",
    siteIndex: 1,
    siteIndex2: 1,
    hubbardVValue: 1.0,
};

export class HubbardVContextProvider extends HubbardUContextProvider {
    get defaultData() {
        const firstElementOrbitals = this.getValenceOrbitals(this.firstSpecies);
        const secondElementOrbitals = this.getValenceOrbitals(this.secondSpecies);
        return [
            {
                ...defaultHubbardConfig,
                atomicSpecies: this.firstSpecies,
                atomicSpecies2: this.secondSpecies,
                siteIndex2:
                    this.uniqueElementsWithLabels?.length > 1 ? 2 : defaultHubbardConfig.siteIndex2,
                atomicOrbital:
                    firstElementOrbitals.length > 0
                        ? firstElementOrbitals[firstElementOrbitals.length - 1]
                        : defaultHubbardConfig.atomicOrbital,
                atomicOrbital2:
                    secondElementOrbitals.length > 0
                        ? secondElementOrbitals[secondElementOrbitals.length - 1]
                        : defaultHubbardConfig.atomicOrbital2,
            },
        ];
    }

    get firstSpecies() {
        return this.firstElement;
    }

    get secondSpecies() {
        return this.uniqueElementsWithLabels?.length > 1
            ? this.uniqueElementsWithLabels[1]
            : this.firstSpecies;
    }

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: true,
                removable: true,
            },
            items: {
                atomicSpecies: this.defaultFieldStyles,
                atomicOrbital: this.defaultFieldStyles,
                atomicSpecies2: this.defaultFieldStyles,
                atomicOrbital2: this.defaultFieldStyles,
                siteIndex: this.defaultFieldStyles,
                siteIndex2: this.defaultFieldStyles,
                hubbardVValue: this.defaultFieldStyles,
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            title: "",
            description: "Hubbard V parameters for DFT+U+V calculation.",
            type: "array",
            items: {
                type: "object",
                properties: {
                    atomicSpecies: {
                        type: "string",
                        title: "Species 1",
                        enum: this.uniqueElementsWithLabels,
                        default: this.firstSpecies,
                    },
                    siteIndex: {
                        type: "integer",
                        title: "Site no 1",
                        default: defaultHubbardConfig.siteIndex,
                    },
                    atomicOrbital: {
                        type: "string",
                        title: "Orbital 1",
                    },
                    atomicSpecies2: {
                        type: "string",
                        title: "Species 2",
                        enum: this.uniqueElementsWithLabels,
                        default: this.secondSpecies,
                    },
                    siteIndex2: {
                        type: "integer",
                        title: "Site no 2",
                        default:
                            this.uniqueElementsWithLabels?.length > 1
                                ? 2
                                : defaultHubbardConfig.siteIndex2,
                    },
                    atomicOrbital2: {
                        type: "string",
                        title: "Orbital 2",
                    },
                    hubbardVValue: {
                        type: "number",
                        title: "V (eV)",
                        default: defaultHubbardConfig.hubbardVValue,
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
                                        enum:
                                            orbitals.length > 0
                                                ? orbitals
                                                : this.orbitalListByStability,
                                        default:
                                            orbitals.length > 0
                                                ? orbitals[orbitals.length - 1]
                                                : defaultHubbardConfig.atomicOrbital,
                                    },
                                },
                            };
                        }),
                    },
                    atomicSpecies2: {
                        oneOf: this.uniqueElementsWithLabels.map((elementWithLabel) => {
                            const orbitals = this.getValenceOrbitals(
                                this.getElementSymbol(elementWithLabel),
                            );
                            return {
                                properties: {
                                    atomicSpecies2: {
                                        enum: [elementWithLabel],
                                    },
                                    atomicOrbital2: {
                                        enum:
                                            orbitals.length > 0
                                                ? orbitals
                                                : this.orbitalListByStability,
                                        default:
                                            orbitals.length > 0
                                                ? orbitals[orbitals.length - 1]
                                                : defaultHubbardConfig.atomicOrbital2,
                                    },
                                },
                            };
                        }),
                    },
                },
            },
            minItems: 1,
        };
    }
}
