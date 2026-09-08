import { type ContextItemSchema } from "@mat3ra/esse/dist/js/types";
import { setupNunjucksEnvironment } from "@mat3ra/standata";
import { Utils } from "@mat3ra/utils";
import nunjucks from "nunjucks";

export type UnitContext = ContextItemSchema[];
export type ContextName = ContextItemSchema["name"];
export type ContextExtraData = ContextItemSchema["extraData"];
export type ContextData = ContextItemSchema["data"];

export type ContextItem<
    D extends ContextData = ContextData,
    ED extends ContextExtraData = ContextExtraData,
> = {
    data?: D;
    extraData?: ED;
    isEdited?: boolean;
};

/** Context item with template-facing `data` (may be a superset of persisted `S["data"]`). */
export type ContextItemForRendering<
    S extends ContextItemSchema = ContextItemSchema,
    R = S["data"],
> = Omit<S, "data"> & { data: R };

export type Domain = "executable" | "important";

export type EntityName = "unit" | "subworkflow";

/** Minimal bound for provider external context; the full contract is ExternalContext in context/providers/index.ts */
export type BaseExternalContext = object;

/** Application-agnostic external context, set by Workflow / Subworkflow for any provider to consume. */
export type JobExternalContext = {
    jobHasParent: boolean;
};

export type WorkflowExternalContext = {
    workflowHasRelaxation: boolean;
};

const nunjucksEnvironment = setupNunjucksEnvironment(new nunjucks.Environment());

const parseNumericStrings = (_key: string, value: unknown) => {
    if (typeof value === "string" && value.trim() !== "" && !Number.isNaN(Number(value))) {
        return Number(value);
    }

    return value;
};

/**
 * Context providers expose three data layers. Keep them separate.
 *
 *   1. Persistent data — user input / defaults. Stored in the DB as the context item's `data`,
 *      and bound to RJSForm as `formData`.
 *        API: `getData()` / `setData(data)` / `getContextItemData()`
 *
 *   2. Patched persistent (rendering) data — persistent data plus fields derived from external
 *      context (material, application, …) that must NOT be stored or shown in the form.
 *      Default = identity, override `patchForRendering` when the template needs extra fields.
 *        API: `getDataForRendering()` / `getContextItemDataForRendering()`
 *      If the rendering shape is a proper superset of the persistent shape, declare it via the
 *      third generic `DataForRendering` and pair with a dedicated ESSE schema (see
 *      `PointsPathFormDataProvider` + `points_path_data_provider_rendering.json`).
 *
 *   3. Full rendering context — `{ ...patchedItemByName, ...externalContext }` fed to Jinja.
 *      Built by `ExecutionUnit.render` from `getContextItemDataForRendering()` results.
 *
 * Rules:
 *   - Never bake derived fields into `setData`/`getDefaultData` just to reach templates.
 *   - `setData` may normalize legacy payloads to the current persistent schema.
 *   - `ExecutionUnit.render` MUST call `getContextItemData()` for `this.context` (DB) and
 *     `getContextItemDataForRendering()` for `this.renderingContext` (Jinja) separately.
 *
 * See `.cursor/rules/context-provider-data-layers.mdc` for the full guidance.
 */
abstract class ContextProvider<
    S extends ContextItemSchema = ContextItemSchema,
    EC extends BaseExternalContext = BaseExternalContext,
    /** Data passed to templates; defaults to persisted `S["data"]` (identity `patchForRendering`). */
    DataForRendering = S["data"],
> {
    abstract name: S["name"];

    abstract readonly domain: Domain;

    abstract readonly entityName: EntityName;

    protected abstract getDefaultData(): S["data"];

    protected data?: S["data"];

    /** Seeded from persisted context; material providers refresh via `updateMaterialHash()`. */
    extraData?: S["extraData"];

    readonly externalContext: EC;

    isEdited: boolean;

    constructor(contextItem: Partial<S>, externalContext: EC) {
        this.externalContext = externalContext;
        this.isEdited = contextItem.isEdited || false;

        if (contextItem.data !== undefined) {
            this.data = Utils.clone.deepClone(contextItem.data);
        }

        if (contextItem.extraData !== undefined) {
            this.extraData = Utils.clone.deepClone(contextItem.extraData);
        }
    }

    setIsEdited(isEdited: boolean) {
        this.isEdited = isEdited;
    }

    getData(): S["data"] {
        if (this.data !== undefined) {
            return Utils.clone.deepClone(this.data);
        }
        return this.getDefaultData();
    }

    setData(data: S["data"]) {
        this.data = Utils.clone.deepClone(data);
    }

    /**
     * Derive template-facing `data` from persisted `data`. Override when the template needs fields
     * that must not be stored (e.g. coordinates from symmetry point names + lattice).
     */
    // eslint-disable-next-line class-methods-use-this
    protected patchForRendering(data: S["data"]): DataForRendering {
        return data as DataForRendering;
    }

    getDataForRendering(): DataForRendering {
        return this.patchForRendering(this.getData());
    }

    getContextItemData(): S {
        return {
            name: this.name,
            isEdited: this.isEdited,
            data: this.getData(),
            extraData: this.extraData,
        } as S;
    }

    getContextItemDataForRendering(): ContextItemForRendering<S, DataForRendering> {
        return {
            ...this.getContextItemData(),
            data: this.getDataForRendering(),
        };
    }

    /**
     * Resolves Jinja placeholders in persisted context `data` from `scope.global`.
     */
    renderContext(scopeGlobal: Record<string, unknown>): boolean {
        const data = this.getData();
        const dataJson = JSON.stringify(data);
        const renderedJson = nunjucksEnvironment.renderString(dataJson, scopeGlobal);
        const resolved = JSON.parse(renderedJson, parseNumericStrings) as S["data"];

        if (JSON.stringify(data) === JSON.stringify(resolved)) {
            return false;
        }

        this.setData(resolved);
        this.setIsEdited(true);

        return true;
    }

    /**
     * Helper method to find a context item from a unit context array by name.
     * Returns a partial schema object that can be safely passed to constructors.
     */
    protected static findContextItem<S extends ContextItemSchema>(
        unitContext: UnitContext,
        contextName: ContextName,
    ): Partial<S> {
        const item = unitContext.find((item): item is S => item.name === contextName);
        return item || {};
    }
}

export default ContextProvider;
