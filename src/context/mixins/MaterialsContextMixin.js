// import type { ContextProvider } from "@mat3ra/code/dist/js/context";
// import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
// import type { OrderedInMemoryEntityInSet } from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
// import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
// import type { MaterialMixin } from "@mat3ra/made/dist/js/materialMixin";

// export type MaterialsContextMixinType = {
//     materials: (MaterialMixin & InMemoryEntity & OrderedInMemoryEntityInSet)[];
//     initMaterialsContextMixin: () => void;
// };

export function materialsContextMixin(
    // item: ContextProvider & { _materials: (MaterialMixin & InMemoryEntity)[] },
    item,
) {
    const properties = {
        get materials() {
            return this._materials;
        },
        initMaterialsContextMixin() {
            // @ts-ignore
            const materials = this.config.context?.materials;
            // @ts-ignore
            if (!this.constructor.Material) {
                throw Error("MaterialsContextMixin: Material is undefined");
            }
            this._materials =
                materials && materials.length
                    ? materials
                    : // @ts-ignore
                      [this.constructor.Material.createDefault()];
        },
    }; // as MaterialsContextMixinType & typeof item;

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
