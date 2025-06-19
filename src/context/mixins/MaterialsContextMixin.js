import { contextProvidersGlobalSettings } from "../providers/settings";

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
                    : [contextProvidersGlobalSettings.Material.createDefault()];
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
