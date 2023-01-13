import { Application } from "@exabyte-io/ade.js";
import { ApplicationContextMixin, ContextProvider } from "@exabyte-io/code.js/dist/context";
import { mix } from "mixwith";

export class MLSettingsContextProvider extends mix(ContextProvider).with(ApplicationContextMixin) {
    static Application = Application;

    // eslint-disable-next-line class-methods-use-this
    get uiSchema() {
        return {
            target_column_name: {},
            problem_category: {},
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return {
            target_column_name: "target",
            problem_category: "regression",
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: " ",
            description: "Settings important to machine learning runs.",
            type: "object",
            properties: {
                target_column_name: {
                    type: "string",
                    default: this.defaultData.target_column_name,
                },
                problem_category: {
                    type: "string",
                    default: this.defaultData.problem_category,
                    enum: ["regression", "classification", "clustering"],
                },
            },
        };
    }
}
