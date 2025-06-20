import CryptoJS from "crypto-js";

export function methodDataContextMixin(item) {
    const properties = {
        _methodData: undefined,

        isEdited: false,

        methodDataHash: undefined,

        extraData: undefined,

        initMethodDataContextMixin() {
            this._methodData = (this.config.context && this.config.context.methodData) || {};
            this.isEdited = Boolean(this.config.isEdited);
        },

        /* @summary Replace the logic in constructor with this in order to enable passing `methodDataHash` between
         *          subsequent initializations of the derived class. Not used at present and kept for the record.
         */
        _initMethodDataHash() {
            this.methodDataHash = CryptoJS.MD5(JSON.stringify(this.methodData)).toString();
            this.extraData = { methodDataHash: this.methodDataHash };
            if (!this._methodData) {
                this._methodData = {};
                this.isEdited = false;
                // Commented out to reduce effect on performance. Uncomment for debugging purposes.
                // TODO: remove on next refactoring or convert to log
                // console.warn("MethodDataContextMixin: methodData is undefined or null");
            } else if (this.isMethodDataUpdated) {
                this.isEdited = false;
            } else {
                // eslint-disable-next-line no-undef
                this.isEdited = config.isEdited;
            }
        },

        get methodData() {
            return this._methodData;
        },

        get isMethodDataUpdated() {
            return Boolean(this.extraData && this.extraData.methodDataHash !== this.methodDataHash);
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
