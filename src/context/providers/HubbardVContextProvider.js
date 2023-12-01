import { HubbardUContextProvider } from "./HubbardUContextProvider";

const defaultHubbardConfig = {
    atomicSpecies: "",
    atomicOrbital: "2p",
    atomicSpecies2: "",
    atomicOrbital2: "2p",
    siteIndex: 1,
    siteIndex2: 1,
    hubbardVValue: 1.0,
};

export class HubbardVContextProvider extends HubbardUContextProvider {
    get defaultData() {
        return [{ ...defaultHubbardConfig, atomicSpecies2: this.uniqueElements }];
    }

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: true,
                removable: true,
            },
            title: {
                "ui:classNames": "col-xs-12",
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
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U+V calculation.",
            type: "array",
            items: {
                type: "object",
                properties: {
                    atomicSpecies: {
                        type: "string",
                        title: "Species 1",
                        enum: this.uniqueElements,
                        default: this.uniqueElements?.length > 0 ? this.uniqueElements[0] : "",
                    },
                    atomicOrbital: {
                        type: "string",
                        title: "Orbital 1",
                        enum: [
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
                        ],
                        default: defaultHubbardConfig.atomicOrbital,
                    },
                    atomicSpecies2: {
                        type: "string",
                        title: "Species 2",
                        enum: this.uniqueElements,
                        default: this.uniqueElements?.length > 0 ? this.uniqueElements[0] : "",
                    },
                    atomicOrbital2: {
                        type: "string",
                        title: "Orbital 2",
                        enum: [
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
                        ],
                        default: defaultHubbardConfig.atomicOrbital,
                    },
                    siteIndex: {
                        type: "number",
                        title: "Site no 1",
                        default: defaultHubbardConfig.siteIndex,
                    },
                    siteIndex2: {
                        type: "number",
                        title: "Site no 2",
                        default: defaultHubbardConfig.siteIndex,
                    },
                    hubbardVValue: {
                        type: "number",
                        title: "V (eV)",
                        default: defaultHubbardConfig.hubbardVValue,
                    },
                },
            },
            minItems: 1,
        };
    }
}
