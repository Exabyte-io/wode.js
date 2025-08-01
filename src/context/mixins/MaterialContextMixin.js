import { globalSettings } from "../providers/settings";

export function materialContextMixin(item) {
    const properties = {
        _material: undefined,

        updateMaterialHash() {
            if (this.isEditedIsSetToFalseOnMaterialUpdate) this.isEdited = false;
            this.extraData = { materialHash: this.material.hash };
        },

        // Workaround: Material.createDefault() used to initiate workflow reducer and hence here too
        //  does not have an id. Here we catch when such material is used and avoid resetting isEdited
        get isMaterialCreatedDefault() {
            return !this.material.id;
        },

        get isMaterialUpdated() {
            return Boolean(this.extraData && this.extraData.materialHash !== this.material.hash);
        },

        get material() {
            if (!this._material) {
                throw new Error("Material is not set");
            }
            return this._material;
        },

        initMaterialContextMixin() {
            this._material = this.config.context && this.config.context.material;

            if (!this._material) {
                this._material = globalSettings.Material.createDefault();
            }

            this.updateMaterialHash();
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
