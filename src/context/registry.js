import { ContextProviderRegistry as clsInstance } from "@exabyte-io/ade.js";
import { extendAndPatchRegistry } from "@exabyte-io/code.js/dist/context";

import { wodeProviders } from "./providers";

export const ContextProviderRegistry = extendAndPatchRegistry(clsInstance, wodeProviders);
