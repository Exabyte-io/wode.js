import { ContextProviderRegistry as clsInstance } from "@exabyte-io/ade.js";
import { extendAndPatchRegistry } from "@mat3ra/code/dist/js/context";

import { wodeProviders } from "./providers";

export const ContextProviderRegistry = extendAndPatchRegistry(clsInstance, wodeProviders);
