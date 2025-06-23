import JSONSchemaFormDataProvider from "@exabyte-io/ade.js/dist/js/context/JSONSchemaFormDataProvider";

import { materialContextMixin } from "../mixins/MaterialContextMixin";

const defaultHubbardConfig = {
    atomicSpecies: "",
    atomicOrbital: "2p",
    hubbardUValue: 1.0,
};

export class HubbardUContextProvider extends JSONSchemaFormDataProvider {
    constructor(config) {
        super(config);

        this.initMaterialContextMixin();

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

materialContextMixin(HubbardUContextProvider.prototype);
