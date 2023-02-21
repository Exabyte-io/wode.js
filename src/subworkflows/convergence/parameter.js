/* eslint-disable class-methods-use-this */
import lodash from "lodash";

export class ConvergenceParameter {
    constructor({ name, initialValue }) {
        this.name = name;
        this._inititalValue = initialValue;
    }

    /**
     * Getter for initial value as string.
     * Note: this will be used in assignment unit.
     * @return {string}
     */
    get initialValue() {
        if (!lodash.isString(this._inititalValue)) {
            return `${this._inititalValue}`;
        }
        return this._inititalValue;
    }

    /**
     * @summary Defines how to increment the parameter.
     * @return {string} - increment operation used in assignment unit
     */
    get increment() {
        return ""; // overwrite in derived class
    }

    /**
     * Defines content for updating the unit context
     * @return {Object}
     */
    get unitContext() {
        return {};
    }

    /**
     * Defines content for updating the subworkflowContext
     * @return {Object}
     */
    get subworkflowContext() {
        return {};
    }

    /**
     * Defines value once convergence is reached (for 'exit' unit).
     * Note: This is used in assignment unit and most often the variable will be assigned to itself.
     * @return {string}
     */
    get finalValue() {
        return `${this.name}`;
    }
}
