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

        /**
         * @summary Find cutoff values corresponding to wavefunction or charge
         * density with accuracy level in the pseudo method data
         * @param data {Object}: the data object from pseudo method data
         * @param cutoffEntity {String}: "wavefunction", or "density"
         * @param accuracyLevel {String}: "standard", "high", or "low"
         * @return {number}: if cutoff value present returns value else return 0
         */
        _cutoffFromPseudoData(data, cutoffEntity, accuracyLevel = "standard") {
            const cutoff = data?.cutoffs?.[cutoffEntity] || [];
            return cutoff.find((obj) => obj.accuracy_level === accuracyLevel)?.value ?? 0;
        },

        /**
         * @summary Set highest cutoff of all elements present, in case wavefunction
         * cutoff is present but no charge density cutoff, set charge density cutoff
         * to 4 or 8 times that of wavefunction cutoff
         * @return {Array<number>}: tuple of wavefunction and density cutoffs
         */
        get highestCutoffsFromPseudos() {
            let ecutwfc = 0;
            let ecutrho = 0;
            const pseudos = this.methodData?.pseudo || [];

            pseudos.forEach((data) => {
                const ecutwfcStandard = this._cutoffFromPseudoData(data, "wavefunction");
                const ecutrhoStandard = this._cutoffFromPseudoData(data, "density");
                // set the highest cutoff of all elements
                ecutwfc = Math.max(ecutwfc, ecutwfcStandard);

                if (ecutrhoStandard > ecutrho) {
                    ecutrho = ecutrhoStandard;
                } else {
                    // if rho cutoff is not present, set it based on wfc cutoff
                    // if it is ultrasoft pseudopotential set rho cutoff 8 times
                    // that of wfc cutoff, otherwise 4 times that of wfc cutoff
                    const rhoMultiplier = data?.type === "us" ? 8 : 4;
                    ecutrho = Math.max(ecutrho, ecutwfc * rhoMultiplier);
                }
            });

            return [ecutwfc, ecutrho];
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
