import { type Taggable } from "@mat3ra/code/dist/js/entity/mixins/TaggableMixin";
import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import type {
    ExecutionUnitInputItemSchema,
    ExecutionUnitSchema,
    TemplateSchema,
} from "@mat3ra/esse/dist/js/types";
import { ApplicationRegistry, applicationVersionSatisfiesSupportedRange } from "@mat3ra/standata";
import { Utils } from "@mat3ra/utils";

import {
    type AnyContextProvider,
    type ExternalContext,
    createProvider,
} from "../context/providers";
import type ConvergenceParameter from "../convergence/ConvergenceParameter";
import { UnitType } from "../enums";
import {
    type ExecutionUnitSchemaMixin,
    executionUnitSchemaMixin,
} from "../generated/ExecutionUnitSchemaMixin";
import BaseUnit from "./BaseUnit";
import ExecutionUnitInput from "./ExecutionUnitInput";

type Schema = ExecutionUnitSchema;

export type ExecutionUnitConfig = Omit<Partial<Schema>, "application"> & SetApplicationProps;

type SetApplicationProps = Pick<Schema, "application"> &
    Pick<Partial<Schema>, "executable" | "flavor"> &
    SetExecutableProps;

type SetExecutableProps = {
    executableName?: string;
    flavorName?: string;
};

interface ExecutionUnit extends ExecutionUnitSchemaMixin, Taggable {}

class ExecutionUnit extends BaseUnit<Schema> implements Schema {
    inputInstances: ExecutionUnitInput[] = [];

    renderingContext: Partial<ExternalContext> = {};

    contextProvidersInstances: AnyContextProvider[] = [];

    static get jsonSchema() {
        return JSONSchemasInterface.getSchemaById("workflow/unit/execution");
    }

    constructor(config: ExecutionUnitConfig) {
        const schema = {
            name: UnitType.execution,
            type: UnitType.execution as Schema["type"],
            input: [],
            results: [],
            preProcessors: [],
            postProcessors: [],
            monitors: [],
            context: [],
            ...config,
        };
        super(schema);

        this.setApplication(config);

        this.name = this.name || this.flavor.name || "";
    }

    setApplication({
        application,
        executable,
        flavor,
        executableName,
        flavorName,
    }: SetApplicationProps) {
        const currentExecutable = this._json.executable;
        const currentFlavor = this._json.flavor;

        this.application = application;
        this.setExecutable({
            executableName: executableName ?? executable?.name ?? currentExecutable?.name,
            flavorName: flavorName ?? flavor?.name ?? currentFlavor?.name,
        });
    }

    setExecutable({ executableName, flavorName }: SetExecutableProps) {
        const executable = new ApplicationRegistry()
            .getExecutablesByApplication(this.application)
            .find((executable) => {
                return executableName ? executable.name === executableName : executable.isDefault;
            });

        if (!executable) {
            throw new Error(`Executable ${executableName} not found`);
        }

        this.executable = executable;
        this.setFlavor(flavorName);
    }

    setFlavor(flavorName?: string) {
        const flavor = new ApplicationRegistry()
            .getFlavorsByApplicationExecutable(this.application, this.executable)
            .find((flavor) => (flavorName ? flavor.name === flavorName : flavor.isDefault));

        if (!flavor) {
            throw new Error(`Flavor ${flavorName} not found`);
        }

        this.defaultResults = flavor.results;
        this.defaultMonitors = flavor.monitors;
        this.defaultPreProcessors = flavor.preProcessors;
        this.defaultPostProcessors = flavor.postProcessors;

        // flavor is missing on the first run, so do not use getter this.flavor with requiredProperty
        const previousFlavor = this._json.flavor;
        const isSameFlavor = previousFlavor?.name === flavor.name;

        if (!isSameFlavor) {
            this.results = flavor.results;
            this.monitors = flavor.monitors;
            this.preProcessors = flavor.preProcessors;
            this.postProcessors = flavor.postProcessors;
        }

        this.flavor = flavor;
        // Persisted input belongs to the flavor it was saved with: another flavor's files are a
        // different set of templates under names that often collide (`script.py`,
        // `requirements.txt`), so keeping them would silently serve the previous flavor's content
        // - and leave slots the new flavor adds (e.g. `utils.py`) to be filled by whatever row sits
        // in the same position. Reuse persisted rows only when the flavor is unchanged, i.e. when
        // reconstructing a unit from saved JSON rather than switching flavor in the UI.
        this.setDefaultInput({ reusePersistedInput: isSameFlavor });
    }

    /**
     * Persisted `input[].template` must match the current application/executable (and optional
     * applicationVersion). Otherwise the stored template is stale, and we take the default from
     * ApplicationRegistry.
     */
    private isPersistedInputItemCompatible(item: ExecutionUnitInputItemSchema): boolean {
        const { template } = item;

        if (
            template.applicationName !== this.application.name ||
            template.executableName !== this.executable.name
        ) {
            return false;
        }

        if (
            !applicationVersionSatisfiesSupportedRange(
                this.application.version,
                template.applicationVersion ?? "",
            )
        ) {
            return false;
        }

        return true;
    }

    /**
     * Pair each input slot of the current flavor with a row persisted in saved workflow JSON. A slot
     * first claims the persisted row carrying the same `template.name`; slots still unpaired then fall
     * back to the row in the same position, which is how a row the user renamed in the UI is still
     * recognized. Every persisted row is claimed by at most one slot, so a flavor with more input files
     * than the saved JSON has rows cannot end up with the same file twice.
     */
    private matchPersistedInput(
        driverTemplates: TemplateSchema[],
    ): (ExecutionUnitInputItemSchema | undefined)[] {
        const persisted = Array.isArray(this.input) ? this.input : [];
        const claimedIndices = new Set<number>();

        const matchedIndices = driverTemplates.map((driverTemplate) => {
            const index = persisted.findIndex(
                (item, i) => !claimedIndices.has(i) && item?.template?.name === driverTemplate.name,
            );

            if (index !== -1) claimedIndices.add(index);

            return index;
        });

        return matchedIndices.map((matchedIndex, index) => {
            if (matchedIndex !== -1) return persisted[matchedIndex];
            if (claimedIndices.has(index)) return undefined;

            claimedIndices.add(index);

            return persisted[index];
        });
    }

    /**
     * Build `inputInstances` from the current flavor’s defaults (`ApplicationRegistry#getInput(application, flavor)`),
     * merged with persisted `this.input` from saved workflow JSON (see `matchPersistedInput`); incompatible or
     * missing rows use the registry template. `render()` then serializes from these instances into `this.input`,
     * so UI and saved JSON stay aligned when Subworkflow re-serializes units after render. Pass
     * `reusePersistedInput: false` to ignore `this.input` and take every slot from the registry, as `setFlavor`
     * does when the flavor changes and the persisted rows describe the previous flavor's files.
     */
    setDefaultInput({ reusePersistedInput = true }: { reusePersistedInput?: boolean } = {}) {
        const driverTemplates = new ApplicationRegistry().getInput(this.application, this.flavor);
        const persistedItems = reusePersistedInput ? this.matchPersistedInput(driverTemplates) : [];

        this.inputInstances = driverTemplates.map((driverTemplate, index) => {
            const persistedItem = persistedItems[index];

            if (persistedItem && this.isPersistedInputItemCompatible(persistedItem)) {
                return new ExecutionUnitInput(persistedItem);
            }

            return ExecutionUnitInput.createFromTemplate(driverTemplate);
        });

        this.input = this.inputInstances.map((input) => input.toJSON());
    }

    render(externalContext: ExternalContext, convergence?: ConvergenceParameter) {
        this.contextProvidersInstances = this.getContextProvidersInstances(
            externalContext,
            convergence,
        );

        this.saveContext(externalContext);
    }

    private getContextProvidersInstances(
        externalContext: ExternalContext,
        convergence?: ConvergenceParameter,
    ) {
        const uniqueContextProviderNames = [
            ...new Set(
                this.input
                    .map((input) => {
                        return input.template.contextProviders.map((provider) => {
                            return provider.name;
                        });
                    })
                    .flat(),
            ),
        ];

        // TODO: kgrid should be abstracted and selected by user
        const parameterToProviderMap = {
            N_k: "kgrid",
            N_k_nonuniform: "kgrid",
        } as const;

        return uniqueContextProviderNames
            .map((name) => {
                return createProvider(name, this.context, externalContext);
            })
            .map((provider) => {
                if (convergence && provider.name === parameterToProviderMap[convergence.name]) {
                    provider.applyConvergenceParameter(convergence);
                }
                return provider;
            });
    }

    savePersistentContext() {
        const persistentItems = this.contextProvidersInstances.map((p) => p.getContextItemData());
        this.context = persistentItems.filter((c) => c.isEdited);
    }

    renderContext(scopeGlobal: Record<string, unknown>) {
        this.contextProvidersInstances.forEach((provider) => {
            provider.renderContext(scopeGlobal);
        });
    }

    saveRenderingContext(externalContext: ExternalContext) {
        // scopeGlobal resolves provider data only; do not pass it to input Jinja templates.
        // eslint-disable-next-line @typescript-eslint/no-unused-vars -- omitted from Jinja context
        const { scopeGlobal, ...renderingExternalContext } = externalContext;
        const renderingItems = this.contextProvidersInstances.map((p) =>
            p.getContextItemDataForRendering(),
        );
        this.renderingContext = {
            ...Object.fromEntries(renderingItems.map((context) => [context.name, context.data])),
            ...renderingExternalContext,
        };
        this.input = this.inputInstances.map((input) => {
            return input.render(this.renderingContext).toJSON();
        });
    }

    saveContext({ scopeGlobal, ...externalContext }: ExternalContext) {
        if (scopeGlobal) {
            this.renderContext(scopeGlobal);
        }

        this.savePersistentContext();
        this.saveRenderingContext(externalContext);
    }

    getHashObject() {
        const { input, flavor, application, executable } = this.toJSON();

        return {
            ...super.getHashObject(),
            application,
            executable,
            flavor,
            input: Utils.hash.calculateHashFromObject(
                input.map(({ template }) => {
                    return Utils.str.removeEmptyLinesFromString(
                        Utils.str.removeCommentsFromSourceCode(template.content),
                    );
                }),
            ),
        };
    }
}

executionUnitSchemaMixin(ExecutionUnit.prototype);

export default ExecutionUnit;
