import { globalSettings } from "../providers/settings";

export function materialsContextMixin(item) {
    const properties = {
        get materials() {
            return this._materials;
        },
        initMaterialsContextMixin() {
            const materials = this.config.context?.materials;
            this._materials =
                materials && materials.length
                    ? materials
                    : [globalSettings.Material.createDefault()];
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
