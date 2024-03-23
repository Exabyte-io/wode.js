import { JSONSchemaFormDataProvider, MaterialContextMixin } from "@exabyte-io/code.js/dist/context";
import { Made } from "@mat3ra/made";
import { mix } from "mixwith";

const defaultHubbardConfig = {
    atomicSpecies: "",
    atomicOrbital: "2p",
    hubbardUValue: 1.0,
};

export class HubbardUContextProvider extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
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
    }

    get defaultData() {
        return [
            {
                ...defaultHubbardConfig,
                atomicSpecies: this.uniqueElements?.length > 0 ? this.uniqueElements[0] : "",
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
            title: {
                "ui:classNames": "col-xs-12",
            },
            items: {
                atomicSpecies: this.defaultFieldStyles,
                atomicOrbital: this.defaultFieldStyles,
                hubbardUValue: this.defaultFieldStyles,
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
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
                        default:
                            this.uniqueElementsWithLabels?.length > 0
                                ? this.uniqueElementsWithLabels[0]
                                : "",
                    },
                    atomicOrbital: {
                        type: "string",
                        title: "Atomic orbital",
                        enum: this.orbitalList,
                        default: defaultHubbardConfig.atomicOrbital,
                    },
                    hubbardUValue: {
                        type: "number",
                        title: "Hubbard U (eV)",
                        default: defaultHubbardConfig.hubbardUValue,
                    },
                },
            },
        };
    }
}
