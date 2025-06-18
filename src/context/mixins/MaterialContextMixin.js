// import type { ContextProvider } from "@mat3ra/code/dist/js/context";
// import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
// import type { OrderedInMemoryEntityInSet } from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
// import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
// import type { MaterialMixin } from "@mat3ra/made/dist/js/materialMixin";
// import type { ApplicationMixin } from "src/js/applicationMixin";

// import Application from "../../application";

// export type MaterialContextMixinType = {
//     assertMaterial: () => void;
//     isEditedIsSetToFalseOnMaterialUpdate?: boolean;
//     updateMaterialHash: () => void;
//     isMaterialCreatedDefault: boolean;
//     isMaterialUpdated: boolean;
//     material: MaterialMixin & InMemoryEntity & OrderedInMemoryEntityInSet;
//     extraData?: {
//         materialHash: string;
//     };
//     initMaterialContextMixin: () => void;
//     _application: ApplicationMixin;
// };

export function materialContextMixin(
    // item: ContextProvider & { _material?: MaterialMixin & InMemoryEntity },
    item,
) {
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
            // @ts-ignore
            if (!this.constructor.Material) {
                throw Error("MaterialContextMixin: Material is undefined");
            }
            this._material = this.config.context && this.config.context.material;
            if (!this._material) {
                this._material = this.constructor.Material.createDefault();
            }
            this.updateMaterialHash();
        },
    }; //  as MaterialContextMixinType & typeof item

    Object.defineProperties(item, Object.getOwnPropertyDescriptors(properties));
}
