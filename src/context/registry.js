import { ContextProviderRegistry as clsInstance } from "@exabyte-io/ade.js";
import _ from "underscore";

import * as allContextProviders from "./providers";

/** Create registry of context providers
 *
 * @summary: Adds context providers to registry and patches them with a static
 * class variable if classes are provided with `classesToPatch`.
 */
export const createRegistry = (classesToPatch = []) => {
    return _.map(allContextProviders, (instance, name) => {
        classesToPatch.forEach((cls) => (instance[cls.name] = cls));
        return clsInstance.addProvider({
            instance,
            name,
        });
    });
};

export const ContextProviderRegistry = createRegistry();
