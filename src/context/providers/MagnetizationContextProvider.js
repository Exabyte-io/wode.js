import { ContextProvider, ModelContextMixin } from "@exabyte-io/code.js/dist/context";
import { Model } from "@exabyte-io/mode.js";
import { mix } from "mixwith";

export class MagnetizationContextProvider extends mix(ContextProvider).with(ModelContextMixin) {
    static Model = Model;

    // eslint-disable-next-line class-methods-use-this
    get uiSchema() {
        return {};
    }

    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return {
            totalMagnetization: 0.0,
            nSpin: 1,
        };
    }

    // eslint-disable-next-line class-methods-use-this
    transformData(data) {
        return {
            ...data,
            nSpin: data.totalMagnetization > 0.0 ? 2 : 1,
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: " ",
            description: "Magnetization value.",
            type: "object",
            properties: {
                totalMagnetization: {
                    type: "number",
                    default: this.defaultData.totalMagnetization,
                },
            },
        };
    }
}
