import ContextProvider from "@exabyte-io/ade.js/dist/js/context/ContextProvider";

import { applicationContextMixin } from "../mixins/ApplicationContextMixin";

export class MLTrainTestSplitContextProvider extends ContextProvider {
    constructor(config) {
        super(config);
        this.initApplicationContextMixin();
    }

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
            fraction_held_as_test_set: 0.2,
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            title: " ",
            description:
                "Fraction held as the test set. For example, a value of 0.2 corresponds to an 80/20 train/test split.",
            type: "object",
            properties: {
                fraction_held_as_test_set: {
                    type: "number",
                    default: this.defaultData.fraction_held_as_test_set,
                    minimum: 0,
                    maximum: 1,
                },
            },
        };
    }
}

applicationContextMixin(MLTrainTestSplitContextProvider.prototype);
