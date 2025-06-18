export function materialsContextMixin(item) {
    const properties = {
        get materials() {
            return this._materials;
        },
        initMaterialsContextMixin() {
            const materials = this.config.context?.materials;
            if (!this.constructor.Material) {
                throw Error("MaterialsContextMixin: Material is undefined");
            }
            this._materials =
                materials && materials.length
                    ? materials
                    : [this.constructor.Material.createDefault()];
        },
    };

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
