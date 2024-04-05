import { HubbardUContextProvider } from "./HubbardUContextProvider";

const defaultHubbardConfig = {
    hubbardUValue: 1.0,
};

export class HubbardContextProviderLegacy extends HubbardUContextProvider {
    get defaultData() {
        return [
            {
                ...defaultHubbardConfig,
                atomicSpecies: this.firstElement,
                atomicSpeciesIndex: this.uniqueElementsWithLabels?.length > 0 ? 1 : null,
            },
        ];
    }

    speciesIndexFromSpecies = (species) => {
        return this.uniqueElementsWithLabels?.length > 0
            ? this.uniqueElementsWithLabels.indexOf(species) + 1
            : null;
    };

    transformData = (data) => {
        return data.map((row) => ({
            ...row,
            atomicSpeciesIndex: this.speciesIndexFromSpecies(row.atomicSpecies),
        }));
    };

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: false,
                removable: true,
            },
            items: {
                atomicSpecies: this.defaultFieldStyles,
                atomicSpeciesIndex: { ...this.defaultFieldStyles, "ui:readonly": true },
                hubbardUValue: this.defaultFieldStyles,
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U calculation.",
            type: "array",
            uniqueItems: true,
            minItems: 1,
            items: {
                type: "object",
                properties: {
                    atomicSpecies: {
                        type: "string",
                        title: "Atomic species",
                        enum: this.uniqueElementsWithLabels,
                    },
                    atomicSpeciesIndex: {
                        type: "integer",
                        title: "Species index",
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
