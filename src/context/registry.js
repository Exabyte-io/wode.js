import { ContextProviderRegistry as clsInstance } from "@exabyte-io/ade.js";

import * as allContextProviders from "./providers";

_.map(allContextProviders, (instance, name) => clsInstance.addProvider({
    instance,
    name
}));

export const ContextProviderRegistry = clsInstance;

